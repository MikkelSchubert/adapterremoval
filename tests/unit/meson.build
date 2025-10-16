# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>

unit_tests = executable(
    'unit_tests',
    'adapter_database_test.cpp',
    'adapter_detector_test.cpp',
    'adapter_selector_test.cpp',
    'alignment_test.cpp',
    'argparse_test.cpp',
    'barcodes_test.cpp',
    'buffer_test.cpp',
    'counts_test.cpp',
    'debug_test.cpp',
    'errors_test.cpp',
    'fastq_enc_test.cpp',
    'fastq_test.cpp',
    'json_test.cpp',
    'linereader_test.cpp',
    'logging_test.cpp',
    'main_test.cpp',
    'mathutils_test.cpp',
    'read_group_test.cpp',
    'sequence_sets' / 'adapter_set_test.cpp',
    'sequence_sets' / 'common_test.cpp',
    'sequence_sets' / 'sample_set_test.cpp',
    'sequence_sets' / 'sample_test.cpp',
    'sequence_test.cpp',
    'serializer_test.cpp',
    'strutils_test.cpp',
    'testing_test.cpp',
    'testing.cpp',
    'threading_test.cpp',
    'utilities_test.cpp',
    config_inc,
    link_with: [libcore, libsimd],
    dependencies: [libisal, libdeflate, libthreads],
    include_directories: [
        '..' / '..' / 'src',
    ],
)

run_target(
    'run-unit-tests',
    command: [
        unit_tests,
        '--use-colour', 'yes',
    ],
)