// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 Mikkel Schubert <mikkelsch@gmail.com>
#include "deflate.hpp"  // declarations
#include "buffer.hpp"   // for buffer
#include "debug.hpp"    // for AR_REQUIRE
#include "errors.hpp"   // for gzip_error
#include "strutils.hpp" // for stringify
#include <array>        // for array
#include <cstddef>      // for size_t
#include <cstdint>      // for uint8_t
#include <isa-l.h>      // IWYU pragma: keep
#include <libdeflate.h> // for libdeflate_crc32
#include <memory>       // for unique_ptr
#include <sstream>      // for ostringstream

namespace adapterremoval {

//! The default buffer size for compression at level ISAL_COMPRESSION_LEVEL
constexpr size_t ISAL_BUFFER_SIZE = ISAL_DEF_LVL1_SMALL;

////////////////////////////////////////////////////////////////////////////////
// isa-l

buffer
gzip_stream_header_isal()
{
  buffer buf(10);

  struct isal_zstream stream = {};
  struct isal_gzip_header header = {};

  isal_gzip_header_init(&header);
  isal_deflate_init(&stream);
  stream.avail_in = 0;
  stream.flush = NO_FLUSH;
  stream.level = ISAL_COMPRESSION_LEVEL;
  stream.gzip_flag = IGZIP_GZIP_NO_HDR;
  stream.next_out = buf.data();
  stream.avail_out = buf.size();

  const auto ret = isal_write_gzip_header(&stream, &header);
  // GCOVR_EXCL_BR_START
  AR_REQUIRE(ret == 0, "buffer was not large enough for gzip header");
  AR_REQUIRE(buf.size() == stream.total_out, "gzip header has changed");
  // GCOVR_EXCL_BR_STOP

  return buf;
}

buffer
gzip_stream_tail_isal(const uint32_t input_crc32, const size_t input_size)
{
  buffer buf;
  buf.reserve(8);
  buf.append_u32(input_crc32);
  buf.append_u32(input_size);

  return buf;
}

void
deflate_block_isal_1(const buffer& src,
                     buffer& dst,
                     const size_t dst_offset,
                     const bool eof)
{
  AR_REQUIRE(dst_offset < dst.size());

  thread_local auto ptr = std::make_unique<isal_zstream>();
  thread_local buffer level_buf;
  if (level_buf.capacity()) {
    isal_deflate_reset(ptr.get());
  } else {
    isal_deflate_stateless_init(ptr.get());
    ptr->flush = FULL_FLUSH;
    ptr->level = ISAL_COMPRESSION_LEVEL;
    ptr->level_buf_size = ISAL_BUFFER_SIZE;
    level_buf.resize(ISAL_BUFFER_SIZE);
    ptr->level_buf = level_buf.data();
  }

  ptr->end_of_stream = eof;
  ptr->avail_in = src.size();
  ptr->next_in = const_cast<uint8_t*>(src.data());
  ptr->next_out = dst.data() + dst_offset;
  ptr->avail_out = dst.size() - dst_offset;

  const auto ec = isal_deflate_stateless(ptr.get());
  // The easily compressible input should fit in a single output block
  AR_REQUIRE(ec == COMP_OK, "isal_deflate_stateless returned " + stringify(ec));
  // Resize the buffer to the actually used size
  dst.resize(dst_offset + ptr->total_out);
}

////////////////////////////////////////////////////////////////////////////////
// libdeflate

namespace {

/** Wrapper for lazy initialization of libdeflate compressors */
struct compressor_wrapper
{
public:
  compressor_wrapper(unsigned compression_level)
    : m_compression_level(compression_level)
  {
  }

  libdeflate_compressor* get()
  {
    if (!m_ptr) {
      m_ptr.reset(libdeflate_alloc_compressor(m_compression_level));
      AR_REQUIRE(m_ptr); // GCOVR_EXCL_BR_LINE
    }

    return m_ptr.get();
  }

private:
  const unsigned m_compression_level;
  std::unique_ptr<libdeflate_compressor, decltype(&libdeflate_free_compressor)>
    m_ptr{ nullptr, libdeflate_free_compressor };
};

} // namespace

void
deflate_block_libdeflate(const buffer& src,
                         buffer& dst,
                         const size_t dst_offset,
                         const unsigned compression_level)
{
  AR_REQUIRE(dst_offset < dst.size());
  AR_REQUIRE(compression_level <= MAX_COMPRESSION_LEVEL);
  thread_local std::array<compressor_wrapper, MAX_COMPRESSION_LEVEL + 1>
    compressors{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

  auto* compressor = compressors.at(compression_level).get();
  auto output_size = ::libdeflate_deflate_compress(compressor,
                                                   src.data(),
                                                   src.size(),
                                                   dst.data() + dst_offset,
                                                   dst.size() - dst_offset);

  // The easily compressible input should fit in a single output block
  AR_REQUIRE(output_size); // GCOVR_EXCL_BR_LINE
  // Resize the buffer to the actually used size
  dst.resize(dst_offset + output_size);
}

////////////////////////////////////////////////////////////////////////////////
// other

void
deflate_block_uncompressed(const buffer& src,
                           buffer& dst,
                           const size_t dst_offset)
{
  AR_REQUIRE(src.size() <= 0xFFFF);

  const auto final_size = src.size() + 5 + dst_offset;
  AR_REQUIRE(final_size <= dst.size());

  dst.resize(dst_offset);
  dst.append_u8(1); // BFINAL=1, BTYPE=00; see RFC1951
  dst.append_u16(src.size());
  dst.append_u16(~src.size());
  dst.append(src);
  dst.resize(final_size);
}

uint32_t
calculate_crc32(const buffer& buf, uint32_t initial)
{
  // libdeflate_crc32 returns 0 if buffer is NULL
  if (buf.size() == 0) {
    return initial;
  }

  return libdeflate_crc32(initial, buf.data(), buf.size());
}

namespace {

[[noreturn]] void
throw_isal_inflate_error(const inflate_state* state,
                         const int error,
                         const size_t* isize,
                         const std::string_view filename)
{
  AR_REQUIRE(state);
  std::ostringstream os;
  if (!filename.empty()) {
    os << "error while decompressing " << shell_escape(filename) << ": ";
  }

  switch (error) {
    case ISAL_END_INPUT:
      os << "file is truncated";
      break;
    case ISAL_OUT_OVERFLOW:
      os << "content is larger than ISIZE";
      break;
    case ISAL_INCORRECT_CHECKSUM:
      if (state->block_state != ISAL_BLOCK_FINISH) {
        os << "incorrect header checksum";
      } else if (isize) {
        if (*isize == state->total_out) {
          os << "incorrect content checksum";
        } else {
          os << "content is smaller than ISIZE";
        }
      } else {
        // isa-l does not distinguish
        os << "incorrect content checksum or ISIZE";
      }
      break;
    case ISAL_INVALID_BLOCK:
    case ISAL_INVALID_SYMBOL:
    case ISAL_INVALID_LOOKBACK:
      os << "invalid deflate block";
      break;
    default:
      // for most part, the exact error is not useful to the user
      os << "not a valid gzip file";
  }

  os << " (error-code = " << error << ")";
  if (error != ISAL_END_INPUT) {
    os << "; file is likely corrupt";
  }

  throw gzip_error(os.str());
}

} // namespace

void
checked_isal_inflate(inflate_state* state,
                     const bool eof,
                     const std::string_view filename)
{
  AR_REQUIRE(state);

  const auto result = isal_inflate(state);
  if (result != ISAL_DECOMP_OK) {
    throw_isal_inflate_error(state, result, nullptr, filename);
  } else if (eof && state->avail_in == 0 &&
             state->block_state != isal_block_state::ISAL_BLOCK_FINISH) {
    std::ostringstream os;
    if (!filename.empty()) {
      os << "error while decompressing " << shell_escape(filename) << ": ";
    }

    os << "file is truncated";

    throw gzip_error(os.str());
  }
}

void
checked_isal_inflate_stateless(inflate_state* state,
                               const size_t isize,
                               const std::string_view filename)
{
  AR_REQUIRE(state);

  const auto result = isal_inflate_stateless(state);
  if (result != ISAL_DECOMP_OK) {
    throw_isal_inflate_error(state, result, &isize, filename);
  }
}

} // namespace adapterremoval
