// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2022 Mikkel Schubert <mikkelsch@gmail.com>
#include "simd.hpp"  // declarations
#include "debug.hpp" // for AR_FAIL

#include "config-ar3.hpp" // for HAVE_SSE2, HAVE_AVX2, HAVE_AVX512...

namespace adapterremoval {

namespace simd {

bool
compare_subsequences_std(size_t& n_mismatches,
                         size_t& n_ambiguous,
                         const char* seq_1,
                         const char* seq_2,
                         size_t length,
                         size_t max_penalty);

bool
compare_subsequences_sse2(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1,
                          const char* seq_2,
                          size_t length,
                          size_t max_penalty);

bool
compare_subsequences_avx2(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1,
                          const char* seq_2,
                          size_t length,
                          size_t max_penalty);

bool
compare_subsequences_avx512(size_t& n_mismatches,
                            size_t& n_ambiguous,
                            const char* seq_1,
                            const char* seq_2,
                            size_t length,
                            size_t max_penalty);

bool
compare_subsequences_neon(size_t& n_mismatches,
                          size_t& n_ambiguous,
                          const char* seq_1,
                          const char* seq_2,
                          size_t length,
                          size_t max_penalty);

std::vector<instruction_set>
supported()
{
  std::vector<instruction_set> choices = { instruction_set::none };

#if HAVE_SSE2
  if (__builtin_cpu_supports("sse2")) {
    choices.push_back(instruction_set::sse2);
  }
#endif

#if HAVE_AVX2
  if (__builtin_cpu_supports("avx2")) {
    choices.push_back(instruction_set::avx2);
  }
#endif

#if HAVE_AVX512
  if (__builtin_cpu_supports("avx512bw")) {
    choices.push_back(instruction_set::avx512);
  }
#endif

#if HAVE_NEON
  choices.push_back(instruction_set::neon);
#endif

  return choices;
}

const char*
name(instruction_set value)
{
  switch (value) {
    case instruction_set::none:
      return "none";
    case instruction_set::sse2:
      return "SSE2";
    case instruction_set::avx2:
      return "AVX2";
    case instruction_set::avx512:
      return "AVX512";
    case instruction_set::neon:
      return "NEON";
    default:
      AR_FAIL("SIMD function not implemented!");
  }
}

size_t
padding(instruction_set value)
{
  switch (value) {
    case instruction_set::none:
      return 0;
    case instruction_set::sse2:
#if HAVE_SSE2
      return 15;
#else
      AR_FAIL("SSE2 not supported!");
#endif
    case instruction_set::avx2:
#if HAVE_AVX2
      return 31;
#else
      AR_FAIL("AVX2 not supported!");
#endif
    case instruction_set::avx512:
#if HAVE_AVX512
      return 0;
#else
      AR_FAIL("AVX512 not supported!");
#endif
    case instruction_set::neon:
#if HAVE_NEON
      return 15;
#else
      AR_FAIL("NEON not supported!");
#endif

    default:
      AR_FAIL("SIMD function not implemented!");
  }
}

compare_subsequences_func
get_compare_subsequences_func(instruction_set is)
{
  switch (is) {
    case instruction_set::none:
      return &compare_subsequences_std;
    case instruction_set::sse2:
#if HAVE_SSE2
      return &compare_subsequences_sse2;
#else
      AR_FAIL("SSE2 not supported!");
#endif
    case instruction_set::avx2:
#if HAVE_AVX2
      return &compare_subsequences_avx2;
#else
      AR_FAIL("AVX2 not supported!");
#endif
    case instruction_set::avx512:
#if HAVE_AVX512
      return &compare_subsequences_avx512;
#else
      AR_FAIL("AVX512 not supported!");
#endif
    case instruction_set::neon:
#if HAVE_NEON
      return &compare_subsequences_neon;
#else
      AR_FAIL("NEON not supported!");
#endif
    default:
      AR_FAIL("SIMD function not implemented!");
  }
}

} // namespace simd

} // namespace adapterremoval
