# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>

project(
    'adapterremoval3',
    'cpp',
    version: '3.0.0-alpha3',
    meson_version: '>=1.2.0',
    default_options: {
        # Disable costly (libc++) assertions in release mode
        'b_ndebug': 'true',
        'cpp_std': 'c++17',
        'debug': 'false',
        'optimization': '2',
        'warning_level': '3',
    },
)

compiler = meson.get_compiler(
    'cpp',
)
compiler_id = compiler.get_id()
compiler_version = compiler.version()

# '' is defined to always be not found
not_found = dependency('', required: false)

################################################################################

add_project_arguments(
    compiler.get_supported_arguments(
        '-fno-rtti',
        '-pedantic-errors',
        '-Wfatal-errors', # Terminate on the first error to simplify debugging
        # Convert certain warnings to errors
        '-Werror=return-type', # Treat missing/wrong return statements as errors
        '-Werror=switch-enum', # Require all enum values in a switch, to catch unhandled cases
        # Additional warnings
        '-Wcast-align',
        '-Wcast-qual',
        '-Wctor-dtor-privacy',
        '-Wdate-time',
        '-Wdeprecated',
        '-Wdisabled-optimization',
        '-Weffc++',
        '-Wformat=2',
        '-Winit-self',
        '-Wold-style-cast',
        '-Woverloaded-virtual',
        '-Wpedantic',
        '-Wredundant-decls',
        '-Wsign-promo',
        '-Wstrict-overflow=2',
        '-Wswitch-default', # Require a default: statement in a switch, to catch corruption
        '-Wundef', # Undefined variables used in #if pragmas
        checked: 'off',
    ),
    language: 'cpp',
)

# Build static binary, forcing linking against static libraries
static = get_option('static')
if static
    add_project_link_arguments(
        compiler.get_supported_link_arguments(
            '-static',
        ),
        language: 'cpp',
    )
endif

# Build static binary, forcing linking against static libraries
libmimalloc = not_found
if get_option('mimalloc')
    libmimalloc = dependency('libmimalloc', 'mimalloc', version: ['>=2', '<3'], static: static)
endif

# Hardening flags based on https://wiki.debian.org/Hardening
harden = get_option('harden')
if harden
    add_project_arguments(
        compiler.get_supported_arguments(
            '-D_FORTIFY_SOURCE=2',
            '-fPIE',
            '-fstack-protector-all',
            '-fstack-clash-protection',
        ),
        language: 'cpp',
    )

    add_project_link_arguments(
        compiler.get_supported_link_arguments(
            '-pie',
            '-Wl,-z,relro',
            '-Wl,-z,now',
            '-Wl,-z,noexecstack',
        ),
        language: 'cpp',
    )
endif

libisal = dependency('libisal', 'isal', version: ['>=2.30.0', '<3'], static: static)
libdeflate = dependency('libdeflate', 'deflate', version: ['>=1', '<2'], static: static)
libthreads = dependency('threads', static: static)

################################################################################

subdir('src')
subdir('docs')
subdir('examples')
subdir('tests')

install_data(
    'README.md',
    install_dir: get_option('datadir') / meson.project_name(),
)

################################################################################
# Helper targets for building static/distro independent binaries

# Name of podman container for static compilation
container_name = 'ar3static'

container_runner = find_program('podman', 'docker', required: false)
if container_runner.found()
    run_target(
        'static-container',
        command: [
            container_runner,
            'build',
            '-f', meson.current_source_dir() / 'Containerfile',
            '-t', container_name,
            meson.current_source_dir(),
        ],
    )

    run_target(
        'static',
        command: [
            container_runner,
            'run',
            '--rm',
            '-t',
            '--mount', 'type=bind,src=@0@/,dst=/host/src/'.format(meson.current_source_dir()),
            '--mount', 'type=bind,src=@0@/,dst=/host/out/'.format(meson.current_build_dir()),
            '--entrypoint', '/usr/bin/make',
            container_name,
            '-C', '/host/src',
            'BUILDDIR=/host/out/static/build',
            'DESTDIR=/host/out/static/install',
            'DEBUG=@0@'.format(get_option('debug')),
            'COVERAGE=@0@'.format(get_option('b_coverage')),
            # Compilation with sanitize flags fails with alpine 3.20.1, but support is
            # left in to avoid giving the false impression that they were enabled
            'SANITIZE=@0@'.format(get_option('b_sanitize') != 'none'),
            'HARDEN=@0@'.format(get_option('harden')),
            'STATIC=true',
            'MIMALLOC=true',
            'setup',
            'tests',
            'install',
        ],
    )
endif