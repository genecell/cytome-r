// cytome: the typed arrays inside a decompressed blob.
//
// A blob holds a numpy array's bytes (`.tobytes()`): the values in the dtype
// the row's `dtype` column names, little-endian unless the name says '>', and
// index arrays as int32 -- or int64, which some writers produced (scipy keeps
// int64 indices through slicing) and which is recognised by the blob's length,
// as the Python reader does. Every array is checked against the sizes the
// chunk records before it is used: a blob of the wrong length or an index out
// of range is an error, never a matrix with values in the wrong place.
#include <Rcpp.h>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
using namespace Rcpp;

typedef std::vector<unsigned char> ByteBuf;
void cytome_inflate(const unsigned char* src, size_t n, std::string method, ByteBuf& out);

static bool host_is_little_endian() {
  const uint16_t one = 1;
  unsigned char b;
  std::memcpy(&b, &one, 1);
  return b == 1;
}

struct DType { char kind; int size; bool big; };

// numpy's names (str(np.dtype(...))) and type codes, with an optional byte-order mark.
static DType parse_dtype(const std::string& name) {
  std::string s = name;
  bool big = false;
  if (!s.empty() && (s[0] == '<' || s[0] == '=' || s[0] == '|')) s = s.substr(1);
  else if (!s.empty() && s[0] == '>') { big = true; s = s.substr(1); }
  static const struct { const char* name; char kind; int size; } table[] = {
    {"bool", 'b', 1}, {"b1", 'b', 1}, {"?", 'b', 1},
    {"int8", 'i', 1}, {"i1", 'i', 1}, {"int16", 'i', 2}, {"i2", 'i', 2},
    {"int32", 'i', 4}, {"i4", 'i', 4}, {"int64", 'i', 8}, {"i8", 'i', 8},
    {"uint8", 'u', 1}, {"u1", 'u', 1}, {"uint16", 'u', 2}, {"u2", 'u', 2},
    {"uint32", 'u', 4}, {"u4", 'u', 4}, {"uint64", 'u', 8}, {"u8", 'u', 8},
    {"float16", 'f', 2}, {"f2", 'f', 2}, {"float32", 'f', 4}, {"f4", 'f', 4},
    {"float64", 'f', 8}, {"f8", 'f', 8}};
  for (size_t k = 0; k < sizeof(table) / sizeof(table[0]); ++k)
    if (s == table[k].name) return DType{table[k].kind, table[k].size, big};
  stop("cytome: unsupported dtype '%s'", name.c_str());
}

// IEEE 754 half precision to double.
static double half_to_double(uint16_t h) {
  int sign = (h >> 15) & 1, expo = (h >> 10) & 0x1F, frac = h & 0x3FF;
  double v;
  if (expo == 0) v = std::ldexp(static_cast<double>(frac), -24);
  else if (expo == 31) v = frac ? NAN : INFINITY;
  else v = std::ldexp(static_cast<double>(frac | 0x400), expo - 25);
  return sign ? -v : v;
}

template <typename T>
static void decode_as(const unsigned char* p, R_xlen_t n, bool swap, double* out) {
  const int size = sizeof(T);
  unsigned char tmp[8];
  for (R_xlen_t k = 0; k < n; ++k, p += size) {
    if (swap) for (int b = 0; b < size; ++b) tmp[b] = p[size - 1 - b];
    else std::memcpy(tmp, p, size);
    T v;
    std::memcpy(&v, tmp, size);
    out[k] = static_cast<double>(v);
  }
}

static const double EXACT = 9007199254740992.0;   // 2^53: past it a double is not exact

static void decode_values_into(const ByteBuf& b, const DType& dt, double* out, R_xlen_t n,
                               const char* what) {
  const bool swap = dt.size > 1 && (dt.big == host_is_little_endian());
  const unsigned char* p = b.data();
  switch (dt.kind) {
  case 'b':
    for (R_xlen_t k = 0; k < n; ++k) out[k] = p[k] != 0;
    return;
  case 'i':
    if (dt.size == 1) return decode_as<int8_t>(p, n, swap, out);
    if (dt.size == 2) return decode_as<int16_t>(p, n, swap, out);
    if (dt.size == 4) return decode_as<int32_t>(p, n, swap, out);
    decode_as<int64_t>(p, n, swap, out);
    break;
  case 'u':
    if (dt.size == 1) return decode_as<uint8_t>(p, n, swap, out);
    if (dt.size == 2) return decode_as<uint16_t>(p, n, swap, out);
    if (dt.size == 4) return decode_as<uint32_t>(p, n, swap, out);
    decode_as<uint64_t>(p, n, swap, out);
    break;
  case 'f':
    if (dt.size == 4) return decode_as<float>(p, n, swap, out);
    if (dt.size == 8) return decode_as<double>(p, n, swap, out);
    for (R_xlen_t k = 0; k < n; ++k, p += 2) {
      uint16_t h = swap ? static_cast<uint16_t>((p[0] << 8) | p[1]) : static_cast<uint16_t>(p[0] | (p[1] << 8));
      out[k] = half_to_double(h);
    }
    return;
  }
  // 64-bit integers: R holds them as doubles, exact only up to 2^53.
  for (R_xlen_t k = 0; k < n; ++k)
    if (std::fabs(out[k]) > EXACT)
      stop("%s: a 64-bit value (%.0f) is past 2^53 and cannot be held exactly in R", what, out[k]);
}

static void inflate_blob(RawVector blob, const std::string& compression, ByteBuf& out) {
  cytome_inflate(reinterpret_cast<const unsigned char*>(blob.begin()),
                 static_cast<size_t>(blob.size()), compression, out);
}

// `n` indices from `b`: int32, or int64 when the blob is twice that long. Every
// index must lie in [0, upper).
static void decode_index_into(const ByteBuf& b, R_xlen_t n, double upper, int* out, const char* what) {
  const size_t bytes = b.size();
  const bool swap = !host_is_little_endian();
  const unsigned char* p = b.data();
  if (bytes == static_cast<size_t>(n) * 4) {
    for (R_xlen_t k = 0; k < n; ++k, p += 4) {
      unsigned char t[4];
      if (swap) { t[0] = p[3]; t[1] = p[2]; t[2] = p[1]; t[3] = p[0]; } else std::memcpy(t, p, 4);
      int32_t v;
      std::memcpy(&v, t, 4);
      if (v < 0 || static_cast<double>(v) >= upper)
        stop("%s: index %d is out of range [0, %.0f)", what, v, upper);
      out[k] = v;
    }
    return;
  }
  if (bytes == static_cast<size_t>(n) * 8 && n > 0) {
    for (R_xlen_t k = 0; k < n; ++k, p += 8) {
      unsigned char t[8];
      if (swap) for (int q = 0; q < 8; ++q) t[q] = p[7 - q]; else std::memcpy(t, p, 8);
      int64_t v;
      std::memcpy(&v, t, 8);
      if (v < 0 || static_cast<double>(v) >= upper || v > INT_MAX)
        stop("%s: index %.0f is out of range [0, %.0f)", what, static_cast<double>(v), upper);
      out[k] = static_cast<int>(v);
    }
    return;
  }
  stop("%s: %.0f bytes for %.0f indices, neither int32 nor int64", what,
       static_cast<double>(bytes), static_cast<double>(n));
}

//' Values of a blob, as doubles
//' @noRd
// [[Rcpp::export]]
NumericVector cytome_decode_values(RawVector blob, std::string compression, std::string dtype,
                                   double n = -1, std::string what = "blob") {
  DType dt = parse_dtype(dtype);
  ByteBuf b;
  inflate_blob(blob, compression, b);
  if (b.size() % dt.size)
    stop("%s: %.0f bytes is not a whole number of %s values", what.c_str(),
         static_cast<double>(b.size()), dtype.c_str());
  R_xlen_t len = static_cast<R_xlen_t>(b.size() / dt.size);
  if (n >= 0 && static_cast<double>(len) != n)
    stop("%s: %.0f %s values where %.0f are recorded", what.c_str(), static_cast<double>(len),
         dtype.c_str(), n);
  NumericVector out(len);
  decode_values_into(b, dt, out.begin(), len, what.c_str());
  return out;
}

//' Indices of a blob (int32 or int64 by its length), each in [0, upper)
//' @noRd
// [[Rcpp::export]]
IntegerVector cytome_decode_index(RawVector blob, std::string compression, double n, double upper,
                                  std::string what = "index blob") {
  ByteBuf b;
  inflate_blob(blob, compression, b);
  IntegerVector out(static_cast<R_xlen_t>(n));
  decode_index_into(b, static_cast<R_xlen_t>(n), upper, out.begin(), what.c_str());
  return out;
}

//' One compressed sparse chunk: `n_major` rows (CSR) or columns (CSC) whose
//' entries index `n_minor` columns or rows. Returns the chunk-local pointer
//' `p` (length n_major + 1, from 0), the 0-based minor indices `i` and the
//' values `x`, after checking that the three arrays agree.
//' @noRd
// [[Rcpp::export]]
List cytome_decode_sparse(RawVector data, RawVector indices, RawVector indptr,
                          std::string compression, std::string dtype, int n_major, double n_minor,
                          std::string what = "chunk") {
  if (n_major < 0) stop("%s: negative extent", what.c_str());
  DType dt = parse_dtype(dtype);
  ByteBuf b;
  inflate_blob(indptr, compression, b);
  IntegerVector p(static_cast<R_xlen_t>(n_major) + 1);
  std::string w = what + " indptr";
  decode_index_into(b, static_cast<R_xlen_t>(n_major) + 1, static_cast<double>(INT_MAX) + 1.0,
                    p.begin(), w.c_str());
  if (p[0] != 0) stop("%s: its index pointer starts at %d, not 0", what.c_str(), p[0]);
  for (int k = 0; k < n_major; ++k)
    if (p[k + 1] < p[k]) stop("%s: its index pointer decreases at %d", what.c_str(), k);
  R_xlen_t nnz = p[n_major];
  inflate_blob(indices, compression, b);
  IntegerVector ii(nnz);
  w = what + " indices";
  decode_index_into(b, nnz, n_minor, ii.begin(), w.c_str());
  inflate_blob(data, compression, b);
  if (b.size() != static_cast<size_t>(nnz) * dt.size)
    stop("%s: %.0f bytes of data for %.0f %s values", what.c_str(), static_cast<double>(b.size()),
         static_cast<double>(nnz), dtype.c_str());
  NumericVector x(nnz);
  decode_values_into(b, dt, x.begin(), nnz, what.c_str());
  return List::create(_["p"] = p, _["i"] = ii, _["x"] = x);
}

//' One fragment chunk: starts (encoding 1: deltas from the previous start, the
//' first from 0), ends (encoding 1: lengths) and cells, all int32. The cells
//' are checked against the cells table by the caller.
//' @noRd
// [[Rcpp::export]]
List cytome_decode_fragments(RawVector starts, RawVector ends, RawVector cells,
                             std::string compression, int encoding,
                             std::string what = "fragment chunk") {
  if (encoding != 0 && encoding != 1)
    stop("%s: unknown fragment encoding %d", what.c_str(), encoding);
  ByteBuf b;
  const bool swap = !host_is_little_endian();
  auto int32s = [&](RawVector blob, std::vector<int32_t>& out, const char* name) {
    inflate_blob(blob, compression, b);
    if (b.size() % 4) stop("%s: %s is not a whole number of int32 values", what.c_str(), name);
    out.resize(b.size() / 4);
    for (size_t k = 0; k < out.size(); ++k) {
      unsigned char t[4];
      const unsigned char* q = b.data() + 4 * k;
      if (swap) { t[0] = q[3]; t[1] = q[2]; t[2] = q[1]; t[3] = q[0]; } else std::memcpy(t, q, 4);
      std::memcpy(&out[k], t, 4);
    }
  };
  std::vector<int32_t> s, e, c;
  int32s(starts, s, "starts");
  int32s(ends, e, "ends");
  int32s(cells, c, "cells");
  if (s.size() != e.size() || s.size() != c.size())
    stop("%s: %.0f starts, %.0f ends and %.0f cells", what.c_str(), static_cast<double>(s.size()),
         static_cast<double>(e.size()), static_cast<double>(c.size()));
  const R_xlen_t n = static_cast<R_xlen_t>(s.size());
  IntegerVector start(n), end(n), cell(n);
  int64_t at = 0;
  for (R_xlen_t k = 0; k < n; ++k) {
    int64_t st = encoding == 1 ? (at += s[k]) : s[k];
    int64_t en = encoding == 1 ? st + e[k] : e[k];
    if (st < 0 || st > INT_MAX || en < st || en > INT_MAX)
      stop("%s: fragment %.0f has coordinates %.0f-%.0f", what.c_str(), static_cast<double>(k),
           static_cast<double>(st), static_cast<double>(en));
    start[k] = static_cast<int>(st);
    end[k] = static_cast<int>(en);
    cell[k] = c[k];
  }
  return List::create(_["start"] = start, _["end"] = end, _["cell_idx"] = cell);
}
