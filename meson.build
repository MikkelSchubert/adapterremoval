project(
    'adapterremoval3',
    'cpp',
    version: '3.0.0-alpha3',
    default_options: {
        'buildtype': 'debugoptimized',
        'cpp_std': 'c++17',
        'optimization': '2',
        'warning_level': '3',
        # Disable costly (libc++) assertions in release mode
        'b_ndebug': 'if-release',
    },
)

compiler = meson.get_compiler(
    'cpp',
)
compiler_id = compiler.get_id()
compiler_version = compiler.version()

################################################################################

add_project_arguments(
    compiler.get_supported_arguments(
        '-fno-rtti',
        '-pedantic-errors',
        '-Wfatal-errors', # Terminate on the first error to simplify debugging
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
        '-Wswitch-enum', # Require all enum values in a switch, to catch unhandled cases
        '-Wundef', # Undefined variables used in #if pragmas
        checked: 'off',
    ),
    # TODO: Remove once Makefile has been deleted
    '-DMESON',
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

libisal = dependency('libisal', version: ['>=2.30.0', '<3'], static: static)
libdeflate = dependency('libdeflate', version: ['>=1', '<2'], static: static)
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
            '--mount', 'type=bind,src=@0@/,dst=/root/adapterremoval/'.format(meson.current_source_dir()),
            '--mount', 'type=bind,src=@0@/,dst=/root/build'.format(meson.current_build_dir()),
            container_name,
            'regression',
            'examples',
            'test',
            'install',
            'BUILD_DIR=/root/build',
            'DESTDIR=/root/build/install',
            'VERBOSE=yes',
            'COLOR=no',
            'STATIC=yes',
        ],
    )
endif

################################################################################

summary('type', get_option('buildtype'), section: 'Build')
summary('static', static, section: 'Build')
summary('harden', harden, section: 'Build')