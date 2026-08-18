// cytome: decompress the per-chunk blobs of a .cytome matrix.
//
// Cytome stores each CSR row-block's data/indices/indptr arrays as a separately
// compressed blob (column `compression` in matrix_chunks). The codecs match
// cytome/io/compression.py exactly:
//   - "zstd"  : a zstd frame (content size in the frame header).
//   - "lz4"   : python lz4.block.compress(store_size=True) => a 4-byte little-endian
//               uncompressed-size header followed by a raw LZ4 block.
//   - "zlib"  : zlib.compress() => zlib-wrapped deflate, no stored size (grow buffer).
#include <Rcpp.h>
#include <lz4.h>
#include <zstd.h>
#include <zlib.h>
#include <cstdint>
#include <cstring>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
RawVector cytome_decompress(RawVector blob, std::string method) {
  const char* src = reinterpret_cast<const char*>(blob.begin());
  size_t srclen = static_cast<size_t>(blob.size());

  if (method == "zstd") {
    unsigned long long dsize = ZSTD_getFrameContentSize(src, srclen);
    if (dsize == ZSTD_CONTENTSIZE_ERROR || dsize == ZSTD_CONTENTSIZE_UNKNOWN)
      stop("cytome: zstd frame content size unknown");
    RawVector out(static_cast<R_xlen_t>(dsize));
    size_t r = ZSTD_decompress(out.begin(), static_cast<size_t>(dsize), src, srclen);
    if (ZSTD_isError(r)) stop("cytome: zstd decompress error: %s", ZSTD_getErrorName(r));
    return out;
  }

  if (method == "lz4") {
    if (srclen < 4) stop("cytome: lz4 blob too short for store_size header");
    uint32_t dsize = 0;
    std::memcpy(&dsize, src, 4);                 // little-endian (cytome runs on x86_64)
    RawVector out(static_cast<R_xlen_t>(dsize));
    int r = LZ4_decompress_safe(src + 4, reinterpret_cast<char*>(out.begin()),
                                static_cast<int>(srclen - 4), static_cast<int>(dsize));
    if (r < 0) stop("cytome: lz4 decompress error (code %d)", r);
    if (static_cast<uint32_t>(r) != dsize) stop("cytome: lz4 size mismatch");
    return out;
  }

  if (method == "zlib") {
    uLongf cap = static_cast<uLongf>(srclen) * 4 + 1024;
    std::vector<Bytef> buf(cap);
    for (int attempt = 0; attempt < 24; ++attempt) {
      uLongf out_len = cap;
      int rc = uncompress(buf.data(), &out_len,
                          reinterpret_cast<const Bytef*>(src), static_cast<uLong>(srclen));
      if (rc == Z_OK) {
        RawVector out(static_cast<R_xlen_t>(out_len));
        std::memcpy(out.begin(), buf.data(), out_len);
        return out;
      }
      if (rc == Z_BUF_ERROR) { cap *= 2; buf.assign(cap, 0); continue; }
      stop("cytome: zlib decompress error (code %d)", rc);
    }
    stop("cytome: zlib output exceeded growth limit");
  }

  stop("cytome: unknown compression method '%s'", method.c_str());
}

// cytome: COMPRESS a blob to match cytome/io/compression.py exactly (inverse of
// cytome_decompress). Used by the native R writer (write_cytome). zstd = frame;
// lz4 = 4-byte little-endian uncompressed-size header + raw LZ4 block (lz4.block
// store_size=True); zlib = zlib.compress().
// [[Rcpp::export]]
RawVector cytome_compress(RawVector blob, std::string method) {
  const char* src = reinterpret_cast<const char*>(blob.begin());
  int srclen = static_cast<int>(blob.size());

  if (method == "zstd") {
    size_t bound = ZSTD_compressBound(srclen);
    std::vector<char> buf(bound);
    size_t r = ZSTD_compress(buf.data(), bound, src, srclen, 3);
    if (ZSTD_isError(r)) stop("cytome: zstd compress error: %s", ZSTD_getErrorName(r));
    RawVector out(static_cast<R_xlen_t>(r));
    std::memcpy(out.begin(), buf.data(), r);
    return out;
  }
  if (method == "lz4") {
    int bound = LZ4_compressBound(srclen);
    std::vector<char> buf(bound);
    int r = LZ4_compress_default(src, buf.data(), srclen, bound);
    if (r <= 0) stop("cytome: lz4 compress error (code %d)", r);
    RawVector out(static_cast<R_xlen_t>(r) + 4);
    uint32_t usize = static_cast<uint32_t>(srclen);          // 4-byte LE size header
    out[0] = usize & 0xFF; out[1] = (usize >> 8) & 0xFF;
    out[2] = (usize >> 16) & 0xFF; out[3] = (usize >> 24) & 0xFF;
    std::memcpy(out.begin() + 4, buf.data(), r);
    return out;
  }
  if (method == "zlib") {
    uLongf bound = compressBound(srclen);
    std::vector<Bytef> buf(bound);
    int rc = compress2(buf.data(), &bound, reinterpret_cast<const Bytef*>(src), srclen, 6);
    if (rc != Z_OK) stop("cytome: zlib compress error (code %d)", rc);
    RawVector out(static_cast<R_xlen_t>(bound));
    std::memcpy(out.begin(), buf.data(), bound);
    return out;
  }
  stop("cytome: unknown compression method '%s'", method.c_str());
}
