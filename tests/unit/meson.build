# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>

unit_tests = executable(
    'unit_tests',
    'alignment_test.cpp',
    'argparse_test.cpp',
    'barcodes_test.cpp',
    'buffer_test.cpp',
    'counts_test.cpp',
    'debug_test.cpp',
    'fastq_test.cpp',
    'json_test.cpp',
    'linereader_test.cpp',
    'logging_test.cpp',
    'main_test.cpp',
    'mathutils_test.cpp',
    'sequence_sets_test.cpp',
    'serializer_test.cpp',
    'strutils_test.cpp',
    'utilities_test.cpp',
    config_inc,
    link_with: [libcore, libsimd],
    dependencies: [libisal, libdeflate, libthreads],
    include_directories: [
        '..' / '..' / 'src',
    ],
)

test(
    'unit_tests',
    unit_tests,
    args: [
        '--invisibles',
        '--use-colour', 'yes',
    ],
    suite: ['unit'],
    verbose: true,
)