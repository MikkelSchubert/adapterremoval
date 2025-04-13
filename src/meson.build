# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>

################################################################################
# Convert HTML containing template specifications to C++ classes
html_template = custom_target(
    'html-template',
    input: ['reports_template.html'],
    output: ['@BASENAME@_html.hpp', '@BASENAME@_html.cpp'],
    depend_files: ['debug.hpp'],
    command: [
        find_program(
            meson.project_source_root() / 'scripts' / 'html_template_to_cpp.py',
        ),
        '@INPUT@',
        '@OUTDIR@' / '@BASENAME@_html',
    ],
)

################################################################################

# To be replaced with `simd` module once stable
config = configuration_data(
    {
        'HAVE_SSE2': '0',
        'HAVE_AVX2': '0',
        'HAVE_AVX512': '0',
        'HAVE_NEON': '0',
    },
)

# To be replaced with `simd.check` when support for `avx512bw` is added
if target_machine.cpu_family() == 'x86'
    config.set('HAVE_SSE2', '1')
elif target_machine.cpu_family() == 'x86_64'
    config.set('HAVE_SSE2', '1')
    config.set('HAVE_AVX2', '1')

    # Compiler versions known to support AVX512:
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95483
    # https://github.com/llvm/llvm-project/commit/eae26bf73715994c2bd145f9b6dc3836aa4ffd4f
    if (compiler_id == 'gcc' and compiler_version.version_compare('>=11.0')) or (compiler_id == 'clang' and compiler_version.version_compare('>=8.0'))
        config.set('HAVE_AVX512', '1')
    else
        warning('AVX512 support requires GCC>=11 or clang>=8')
    endif
elif target_machine.cpu_family() == 'aarch64'
    config.set('HAVE_NEON', '1')
endif

libsimd = []
libsimd_config = {
    'HAVE_SSE2': ['sse2', 'simd_sse2.cpp', '-msse2'],
    'HAVE_AVX2': ['avx2', 'simd_avx2.cpp', '-mavx2'],
    'HAVE_AVX512': ['avx512', 'simd_avx512bw.cpp', '-mavx512bw'],
    'HAVE_NEON': ['neon', 'simd_neon.cpp', []],
}

foreach key, simd_config : libsimd_config
    summary(simd_config[0], config.get(key) != '0', section: 'Build')
    if config.get(key) != '0'
        libsimd += [
            static_library(
                simd_config[0],
                simd_config[1],
                cpp_args: simd_config[2],
            ),
        ]
    endif
endforeach

config_inc = configure_file(output: 'config-ar3.hpp', configuration: config)

################################################################################

# Library containing objects shared between the executable and unit tests
libcore = static_library(
    'core',
    'alignment.cpp',
    'argparse.cpp',
    'barcode_table.cpp',
    'debug.cpp',
    'errors.cpp',
    'fastq_enc.cpp',
    'fastq.cpp',
    'json.cpp',
    'linereader.cpp',
    'logging.cpp',
    'managed_io.cpp',
    'mathutils.cpp',
    'read_group.cpp',
    'sequence_sets.cpp',
    'sequence.cpp',
    'serializer.cpp',
    'simd_std.cpp',
    'simd.cpp',
    'strutils.cpp',
    'table_reader.cpp',
    'utilities.cpp',
    config_inc,
    dependencies: [libdeflate, libisal, libmimalloc, libthreads],
)

exe = executable(
    meson.project_name(),
    'adapter_id.cpp',
    'benchmarking.cpp',
    'demultiplexing.cpp',
    'fastq_io.cpp',
    'linereader_joined.cpp',
    'main.cpp',
    'main_adapter_id.cpp',
    'main_adapter_rm.cpp',
    'main_benchmark.cpp',
    'output.cpp',
    'progress.cpp',
    'reports_html.cpp',
    'reports_json.cpp',
    'reports_terminal.cpp',
    'scheduler.cpp',
    'statistics.cpp',
    'timer.cpp',
    'trimming.cpp',
    'userconfig.cpp',
    html_template,
    config_inc,
    link_with: [libcore, libsimd],
    dependencies: [libdeflate, libisal, libmimalloc, libthreads],
    install: true,
)