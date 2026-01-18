# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>

sh = find_program('sh')
sphinx_build = find_program('sphinx-build', required: get_option('docs'))

if sphinx_build.found()
    custom_target(
        'docs',
        input: [
            'conf.py',
            'detailed_overview.rst',
            'examples.rst',
            'getting_started.rst',
            'index.rst',
            'input_and_output.rst',
            'installation.rst',
            'manpage.rst',
            'migrating.rst',
        ],
        output: 'html',
        command: [
            sphinx_build,
            '-M', 'html',
            '@CURRENT_SOURCE_DIR@',
            '@OUTDIR@',
            '-d', '@PRIVATE_DIR@', # do not store doctree files
            '-n', # warn about all missing references
            '-q', # no output on stdout, just warnings on stderr
            '-E', # don't use a saved environment, always read all files
        ],
        install: get_option('docs'),
        install_dir: get_option('datadir') / meson.project_name() / 'docs',
    )

    custom_target(
        'man',
        input: [
            'conf.py',
            'manpage.rst',
        ],
        output: meson.project_name() + '.1',
        # Workaround for meson not accepting subdirectories in output. See
        #   https://github.com/mesonbuild/meson/issues/2320
        command: [
            sh,
            '-c', '@0@ -M man @1@ @PRIVATE_DIR@ -n -q -E -d @PRIVATE_DIR@ && mv -v @2@ @OUTPUT@'.format(
                sphinx_build.full_path(),
                '@CURRENT_SOURCE_DIR@',
                '@PRIVATE_DIR@' / 'man' / (meson.project_name() + '.1'),
            ),
        ],
        install: true,
        install_dir: join_paths(get_option('mandir'), 'man1'),
    )
else
    warning('sphinx-build not found; man-page/documentation cannot be created')
endif