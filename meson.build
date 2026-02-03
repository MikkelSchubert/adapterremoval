# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2024 Mikkel Schubert <mikkelsch@gmail.com>

project(
    'adapterremoval3',
    'cpp',
    version: '3.0.0-alpha3',
    meson_version: '>=0.56.0',
    default_options: [
        # Disable costly (libc++) assertions in release mode
        'b_ndebug=true',
        'cpp_std=c++17',
        'debug=false',
        'optimization=2',
        'warning_level=3',
        'b_lto=true',
        'b_lto_mode=thin',
    ],
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
        # Workaround for systems where _FORTIFY_SOURCE is defined by default
        '-DPROJECT_NAME="@0@"'.format(meson.project_name()),
        '-DPROJECT_VERSION="@0@"'.format(meson.project_version()),
    ),
    language: 'cpp',
)

add_project_arguments(
    compiler.get_supported_arguments(
        # abort when sanitizers trigger
        '-fno-sanitize-recover=all',
        # Treat non-ISO-conforming code and the use of extensions as errors
        '-pedantic-errors',
        # Terminate on the first error to simplify debugging
        '-Wfatal-errors',
        # Convert certain warnings to errors
        '-Werror=return-type', # Treat missing/wrong returns as errors
        '-Werror=switch-enum', # Require all enum values in a switch
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
        '-Wshadow-local',
        '-Wsign-promo',
        '-Wstrict-overflow=2',
        '-Wswitch-default', # Require a default: statement in a switch, to catch corruption
        '-Wundef', # Undefined variables used in #if pragmas
        '-Wnrvo',
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

    section_options = compiler.get_supported_link_arguments(
        '-ffunction-sections',
        '-fdata-sections',
    )

    section_linker_options = compiler.get_supported_link_arguments(
        '-Wl,--gc-sections',
    )

    # -f*sections + --gc-sections reduce the size of the static build by 25%
    if section_options.length() == 2 and section_linker_options.length() == 1
        add_project_arguments(section_options, language: 'cpp')
        add_project_link_arguments(section_linker_options, language: 'cpp')
    endif
endif

# Build static binary, forcing linking against static libraries
libmimalloc = not_found
if get_option('mimalloc')
    libmimalloc = dependency('libmimalloc', version: ['>=2', '<4'], static: static, required: false)
    if not libmimalloc.found()
        libmimalloc = dependency('mimalloc', version: ['>=2', '<4'], static: static)
    endif
endif

# Hardening flags based on https://wiki.debian.org/Hardening
harden = get_option('harden')
if harden
    add_project_arguments(
        compiler.get_supported_arguments(
            # Workaround for systems where _FORTIFY_SOURCE is defined by default
            '-U_FORTIFY_SOURCE',
            '-D_FORTIFY_SOURCE=3',
            '-D_GLIBCXX_ASSERTIONS',
            '-fPIE',
            '-fcf-protection=full',
            '-fstack-clash-protection',
            '-fstack-protector-strong',
            '-ftrivial-auto-var-init=zero',
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

libisal = dependency('libisal', version: ['>=2.30.0', '<3'], static: static, required: false)
if not libisal.found()
    libisal = dependency('isal', version: ['>=2.30.0', '<3'], static: static, required: false)

    # WORKAROUND: Support for older distros lacking cmake/pkgconf support
    if not libisal.found()
        libisal = meson.get_compiler('cpp').find_library('libisal')
    endif
endif

libdeflate = dependency('libdeflate', version: ['>=1', '<2'], static: static, required: false)
if not libdeflate.found()
    libdeflate = dependency('deflate', version: ['>=1', '<2'], static: static, required: false)

    # WORKAROUND: Support for older distros lacking cmake/pkgconf support
    if not libdeflate.found()
        libdeflate = meson.get_compiler('cpp').find_library('libdeflate')
    endif
endif

libthreads = dependency('threads', static: static)

# WORKAROUND: Required for some old distros
libcppfs = meson.get_compiler('cpp').find_library('stdc++fs', required: false)

################################################################################

subdir('src')
subdir('docs')
subdir('examples')
subdir('tests')

install_data(
    'README.md',
    install_dir: get_option('datadir') / meson.project_name(),
)