# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>

sphinx_required = get_option('manpage').enabled() or get_option('docs').enabled()
sphinx_build = find_program('sphinx-build', required: sphinx_required)

if sphinx_build.found()
    custom_target(
        'docs',
        input: [
            'conf.py',
            'detailed_overview.rst',
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
        # Meson >=0.59.0 required for `allowed()`
        install: not get_option('docs').disabled(),
        install_dir: get_option('datadir') / meson.project_name() / 'docs',
    )

    sh = find_program('sh', required: get_option('manpage').enabled())
    if sh.found()
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
            install: not get_option('manpage').disabled(),
            install_dir: join_paths(get_option('mandir'), 'man1'),
        )
    endif
else
    warning('sphinx-build not found; man-page/documentation cannot be created')
endif