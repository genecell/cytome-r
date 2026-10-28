// cytome: the codecs of a .cytome blob.
//
// Each array of a chunk (data, indices, indptr; fragment starts, ends, cells;
// an embedding block) is a separately compressed blob, labelled by the row's
// `compression` column. The codecs match cytome/io/compression.py exactly:
//   - "zstd" : a zstd frame; or zlib / lz4 bytes under that label (see below).
//   - "lz4"  : python lz4.block.compress(store_size=True) => a 4-byte
//              little-endian uncompressed-size header followed by a raw LZ4
//              block. Little-endian on every platform, so it is read byte by
//              byte, never through the host's integer layout.
//   - "zlib" : zlib.compress() => zlib-wrapped deflate, no stored size.
//   - "none" : the bytes themselves (the importer's uncompressed mode).
//
// Sizes: an LZ4 block holds at most LZ4_MAX_INPUT_SIZE bytes (about 2.1 GB),
// a limit of the format that Python's lz4.block shares, so larger blobs are
// refused by name. zstd and zlib take any size R can hold; zlib is driven in
// steps because its length fields are 32-bit on some platforms.
#include <Rcpp.h>
#include <lz4.h>
#include <zstd.h>
#include <zlib.h>
#include <climits>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
using namespace Rcpp;

typedef std::vector<unsigned char> ByteBuf;

// The label is a hint when it says "zstd", as in cytome/io/compression.py:
// without the optional zstandard package the Python writer of cytome 0.3.6
// and earlier compressed with zlib but still recorded "zstd", so the bytes
// decide. "lz4" and "zlib" labels are trusted: lz4's 4-byte size header can
// look like a zlib header.
static bool is_zstd_frame(const unsigned char* s, size_t n) {
  return n >= 4 && s[0] == 0x28 && s[1] == 0xB5 && s[2] == 0x2F && s[3] == 0xFD;
}
static bool is_zlib_stream(const unsigned char* s, size_t n) {
  return n >= 2 && s[0] == 0x78 && (s[1] == 0x01 || s[1] == 0x5E || s[1] == 0x9C || s[1] == 0xDA);
}

static void check_r_size(size_t n) {
  if (n > static_cast<size_t>(R_XLEN_T_MAX))
    stop("cytome: a blob of %.0f bytes is past what R can hold", static_cast<double>(n));
}

// A zstd frame whose header does not record its size (a streaming writer):
// decode it in steps into a growing buffer.
static void zstd_stream(const unsigned char* src, size_t srclen, ByteBuf& out) {
  ZSTD_DStream* ds = ZSTD_createDStream();
  if (!ds) stop("cytome: zstd: cannot allocate a decoder");
  ZSTD_initDStream(ds);
  std::vector<unsigned char> buf(ZSTD_DStreamOutSize());
  ZSTD_inBuffer in = { src, srclen, 0 };
  size_t ret = 1;
  while (in.pos < in.size && ret != 0) {
    ZSTD_outBuffer ob = { buf.data(), buf.size(), 0 };
    ret = ZSTD_decompressStream(ds, &ob, &in);
    if (ZSTD_isError(ret)) {
      ZSTD_freeDStream(ds);
      stop("cytome: zstd decompress error: %s", ZSTD_getErrorName(ret));
    }
    out.insert(out.end(), buf.data(), buf.data() + ob.pos);
  }
  ZSTD_freeDStream(ds);
  if (ret != 0) stop("cytome: zstd frame is truncated");
}

static void zstd_decompress(const unsigned char* src, size_t n, ByteBuf& out) {
  unsigned long long dsize = ZSTD_getFrameContentSize(src, n);
  if (dsize == ZSTD_CONTENTSIZE_ERROR) stop("cytome: a blob labelled zstd is not a zstd frame");
  if (dsize == ZSTD_CONTENTSIZE_UNKNOWN) return zstd_stream(src, n, out);
  check_r_size(static_cast<size_t>(dsize));
  out.resize(static_cast<size_t>(dsize));
  size_t r = ZSTD_decompress(out.empty() ? nullptr : out.data(), out.size(), src, n);
  if (ZSTD_isError(r)) stop("cytome: zstd decompress error: %s", ZSTD_getErrorName(r));
  if (r != out.size()) stop("cytome: zstd frame holds %.0f bytes, its header says %.0f",
                            static_cast<double>(r), static_cast<double>(dsize));
}

static void lz4_decompress(const unsigned char* src, size_t n, ByteBuf& out) {
  if (n < 4) stop("cytome: lz4 blob too short for its size header");
  uint32_t dsize = static_cast<uint32_t>(src[0]) | (static_cast<uint32_t>(src[1]) << 8) |
                   (static_cast<uint32_t>(src[2]) << 16) | (static_cast<uint32_t>(src[3]) << 24);
  if (dsize > static_cast<uint32_t>(LZ4_MAX_INPUT_SIZE))
    stop("cytome: an lz4 blob says it holds %.0f bytes, past the format's limit of %d",
         static_cast<double>(dsize), LZ4_MAX_INPUT_SIZE);
  if (n - 4 > static_cast<size_t>(INT_MAX))
    stop("cytome: an lz4 blob of %.0f bytes is past the format's limit", static_cast<double>(n));
  out.resize(dsize);
  unsigned char empty = 0;
  int r = LZ4_decompress_safe(reinterpret_cast<const char*>(src) + 4,
                              reinterpret_cast<char*>(dsize ? out.data() : &empty),
                              static_cast<int>(n - 4), static_cast<int>(dsize));
  if (r < 0) stop("cytome: lz4 decompress error (code %d)", r);
  if (static_cast<uint32_t>(r) != dsize) stop("cytome: lz4 size mismatch");
}

static void zlib_decompress(const unsigned char* src, size_t n, ByteBuf& out) {
  z_stream zs;
  std::memset(&zs, 0, sizeof(zs));
  if (inflateInit(&zs) != Z_OK) stop("cytome: zlib: cannot allocate a decoder");
  size_t in_pos = 0;
  out.resize(n < 1024 ? 4096 : 4 * n);
  size_t out_pos = 0;
  int rc = Z_OK;
  while (rc != Z_STREAM_END) {
    if (out_pos == out.size()) out.resize(2 * out.size());
    if (zs.avail_in == 0 && in_pos < n) {
      size_t take = n - in_pos;
      if (take > UINT_MAX) take = UINT_MAX;
      zs.next_in = const_cast<Bytef*>(src + in_pos);
      zs.avail_in = static_cast<uInt>(take);
      in_pos += take;
    }
    size_t room = out.size() - out_pos;
    if (room > UINT_MAX) room = UINT_MAX;
    zs.next_out = out.data() + out_pos;
    zs.avail_out = static_cast<uInt>(room);
    rc = inflate(&zs, Z_NO_FLUSH);
    out_pos += room - zs.avail_out;
    if (rc == Z_STREAM_END) break;
    if (rc == Z_BUF_ERROR && zs.avail_in == 0 && in_pos == n && zs.avail_out > 0) {
      inflateEnd(&zs);
      stop("cytome: zlib stream is truncated");
    }
    if (rc != Z_OK && rc != Z_BUF_ERROR) {
      inflateEnd(&zs);
      stop("cytome: zlib decompress error (code %d)", rc);
    }
  }
  // Bytes after the end of the stream are ignored, as Python's zlib.decompress does.
  inflateEnd(&zs);
  check_r_size(out_pos);
  out.resize(out_pos);
}

// Decompress `n` bytes at `src` labelled `method` into `out`. Shared by the
// R-level codec and the typed decoders, so an array is decoded straight from
// the decompressed buffer, without a second copy in R.
void cytome_inflate(const unsigned char* src, size_t n, std::string method, ByteBuf& out) {
  out.clear();
  if (method == "none" || method == "raw" || method.empty()) {
    out.assign(src, src + n);
    return;
  }
  if (method == "zstd" && n >= 4 && !is_zstd_frame(src, n))
    method = is_zlib_stream(src, n) ? "zlib" : "lz4";
  if (method == "zstd") return zstd_decompress(src, n, out);
  if (method == "lz4") return lz4_decompress(src, n, out);
  if (method == "zlib") return zlib_decompress(src, n, out);
  stop("cytome: unknown compression method '%s'", method.c_str());
}

static RawVector to_raw(const ByteBuf& b) {
  RawVector res(static_cast<R_xlen_t>(b.size()));
  if (!b.empty()) std::memcpy(res.begin(), b.data(), b.size());
  return res;
}

// [[Rcpp::export]]
RawVector cytome_decompress(RawVector blob, std::string method) {
  ByteBuf out;
  cytome_inflate(reinterpret_cast<const unsigned char*>(blob.begin()),
                 static_cast<size_t>(blob.size()), method, out);
  return to_raw(out);
}

// cytome: COMPRESS a blob to match cytome/io/compression.py exactly (inverse of
// cytome_decompress). Used by the native R writer (write_cytome). zstd = frame;
// lz4 = 4-byte little-endian uncompressed-size header + raw LZ4 block (lz4.block
// store_size=True); zlib = zlib.compress().
// [[Rcpp::export]]
RawVector cytome_compress(RawVector blob, std::string method) {
  const unsigned char* src = reinterpret_cast<const unsigned char*>(blob.begin());
  size_t n = static_cast<size_t>(blob.size());

  if (method == "zstd") {
    size_t bound = ZSTD_compressBound(n);
    if (ZSTD_isError(bound)) stop("cytome: zstd cannot compress %.0f bytes", static_cast<double>(n));
    ByteBuf buf(bound);
    size_t r = ZSTD_compress(buf.data(), bound, src, n, 3);
    if (ZSTD_isError(r)) stop("cytome: zstd compress error: %s", ZSTD_getErrorName(r));
    buf.resize(r);
    return to_raw(buf);
  }
  if (method == "lz4") {
    if (n > static_cast<size_t>(LZ4_MAX_INPUT_SIZE))
      stop("cytome: an lz4 blob holds at most %d bytes, this one has %.0f; use zstd",
           LZ4_MAX_INPUT_SIZE, static_cast<double>(n));
    int bound = LZ4_compressBound(static_cast<int>(n));
    ByteBuf buf(static_cast<size_t>(bound) + 4);
    unsigned char empty = 0;
    int r = LZ4_compress_default(n ? reinterpret_cast<const char*>(src) : reinterpret_cast<const char*>(&empty),
                                 reinterpret_cast<char*>(buf.data()) + 4, static_cast<int>(n), bound);
    if (r <= 0) stop("cytome: lz4 compress error (code %d)", r);
    uint32_t usize = static_cast<uint32_t>(n);          // 4-byte little-endian size header
    buf[0] = usize & 0xFF; buf[1] = (usize >> 8) & 0xFF;
    buf[2] = (usize >> 16) & 0xFF; buf[3] = (usize >> 24) & 0xFF;
    buf.resize(static_cast<size_t>(r) + 4);
    return to_raw(buf);
  }
  if (method == "zlib") {
    z_stream zs;
    std::memset(&zs, 0, sizeof(zs));
    if (deflateInit(&zs, 6) != Z_OK) stop("cytome: zlib: cannot allocate an encoder");
    ByteBuf out(n / 2 + 1024);
    size_t in_pos = 0, out_pos = 0;
    int rc = Z_OK;
    while (rc != Z_STREAM_END) {
      if (out_pos == out.size()) out.resize(2 * out.size());
      if (zs.avail_in == 0 && in_pos < n) {
        size_t take = n - in_pos;
        if (take > UINT_MAX) take = UINT_MAX;
        zs.next_in = const_cast<Bytef*>(src + in_pos);
        zs.avail_in = static_cast<uInt>(take);
        in_pos += take;
      }
      size_t room = out.size() - out_pos;
      if (room > UINT_MAX) room = UINT_MAX;
      zs.next_out = out.data() + out_pos;
      zs.avail_out = static_cast<uInt>(room);
      rc = deflate(&zs, in_pos == n ? Z_FINISH : Z_NO_FLUSH);   // all input handed over
      out_pos += room - zs.avail_out;
      if (rc != Z_OK && rc != Z_STREAM_END && rc != Z_BUF_ERROR) {
        deflateEnd(&zs);
        stop("cytome: zlib compress error (code %d)", rc);
      }
    }
    deflateEnd(&zs);
    out.resize(out_pos);
    return to_raw(out);
  }
  stop("cytome: unknown compression method '%s'", method.c_str());
}
