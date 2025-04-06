# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>

test_runner_script = meson.project_source_root() / 'scripts' / 'regression_test_runner.py'

regression_test_exe = find_program('uv', required: false, version: '>=0.5.17')
if regression_test_exe.found()
    regression_test_args = ['run', test_runner_script]
else
    regression_test_exe = find_program(test_runner_script)
    regression_test_args = []
endif

custom_target(
    'run-regression-tests',
    command: [regression_test_exe]
    + regression_test_args
    + [
        meson.current_build_dir(),
        meson.current_source_dir(),
        '--use-colors', 'always',
        '--json-schema', meson.project_source_root() / 'schema.json',
        '--executable', exe.full_path(),
        '--source-root', meson.project_source_root(),
    ],
    output: 'regression',
    console: true,
    depends: [exe],
)

custom_target(
    'update-regression-tests',
    command: [regression_test_exe]
    + regression_test_args
    + [
        meson.current_build_dir(),
        meson.current_source_dir(),
        '--create-updated-reference',
        '--use-colors', 'always',
        '--schema-validation-required',
        '--json-schema', meson.project_source_root() / 'schema.json',
        '--executable', exe.full_path(),
        '--source-root', meson.project_source_root(),
    ],
    output: 'update-regression',
    console: true,
    depends: [exe],
)