###############################################################################
# This makefile is provided as a simple convenience

## Optional features; set to 'true' to enable or 'false' to disable:
##   $ make DEBUG=true

# Include coverage instrumentation in build
COVERAGE := false

# Debug build; adds warnings, debugging symbols, and extra STL/POSIX asserts
DEBUG := ${COVERAGE}

# Build and install HTML documentation (https://adapterremoval.readthedocs.org)
DOCS := false

# Enable address and undefined behavior sanitation
SANITIZE := false

# Enable hardening flags
HARDEN := true

# Enable link-time optimizations
LTO := true

# LTO mode; older systems may need to use `make LTO_MODE=default`
LTO_MODE := thin

# Generate statically linked binary
# It is recommended to use the included Containerfile to build the static binary
STATIC := false

# Use/require mimalloc. Mainly intended for static builds, as the musl allocator
# comes with a significant performance cost
MIMALLOC := $(STATIC)

###############################################################################
# (Container for) building static binaries using Alpine

# Podman is preferred (and the fallback if neither is found), since it is more
# likely that the user can execute podman than docker, if both are present
CONTAINER_RUNNER := $(firstword $(shell which podman docker) podman)

CONTAINER_NAME := ar3static

###############################################################################

# Extra meson flags
MESON_OPTIONS :=

# Enable extra checks for STL and POSIX functions
ifeq ($(strip ${DEBUG}), true)
override MESON_OPTIONS += -Db_ndebug=false
endif

ifeq ($(strip ${SANITIZE}), true)
override MESON_OPTIONS += -Db_sanitize=address,undefined
else
ifneq ($(strip ${SANITIZE}), false)
$(error "SANITIZE must be 'true' or 'false', not '${SANITIZE}'")
endif
endif

# Use custom installation prefix instead of the system default
ifneq ($(strip ${PREFIX}), )
override MESON_OPTIONS += --prefix=${PREFIX}
endif

###############################################################################

# Meson commands cannot be run in parallel
.NOTPARALLEL:

# Default meson build directory
BUILDDIR := build
# Location of ninja build-file; used to detected existing setup
NINJAFILE := ${BUILDDIR}/build.ninja

executable: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" adapterremoval3

executables: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" adapterremoval3 unit_tests

clean:
	rm -rf "${BUILDDIR}" ".venv"

clean-coverage:
	test ! -d "${BUILDDIR}" || find "${BUILDDIR}" -name '*.gcda' -delete

coverage: ${NINJAFILE}
	ninja -C "${BUILDDIR}" coverage-text
	cat build/meson-logs/coverage.txt

coverage-xml: ${NINJAFILE}
	ninja -C "${BUILDDIR}" coverage-xml

docs: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" docs man

install: ${NINJAFILE}
	meson install -C "${BUILDDIR}"

regression-tests: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" run-regression-tests

setup ${NINJAFILE}:
	# WORKAROUND: `setup` fails on existing builddirs without `--reconfigure`,
	#             but that option requires existing builddir in v1.0.0 or older
	rm -rf "${BUILDDIR}"
	meson setup "${BUILDDIR}" \
		-Db_coverage=${COVERAGE} \
		-Db_lto=${LTO} \
		-Db_lto_mode=${LTO_MODE} \
		-Ddebug=${DEBUG} \
		-Ddocs=${DOCS} \
		-Dharden=${HARDEN} \
		-Dmimalloc=${MIMALLOC} \
		-Dstatic=${STATIC} \
		${MESON_OPTIONS}

static:
	mkdir -p "${BUILDDIR}"
	# Compilation with sanitize flags fails with alpine, but support is
	# left in to avoid giving the false impression that they were enabled
	"${CONTAINER_RUNNER}" run --rm -t \
		--mount "type=bind,src=${PWD}/,dst=/host/src/" \
		--mount "type=bind,src=${BUILDDIR}/,dst=/host/out/" \
		--entrypoint /usr/bin/make \
		"${CONTAINER_NAME}" \
		-C /host/src \
		BUILDDIR=/host/out/static/build \
		PREFIX=/host/out/static/install \
		DEBUG=${DEBUG} \
		COVERAGE=${COVERAGE} \
		DOCS=${DOCS} \
		SANITIZE=${SANITIZE} \
		HARDEN=${HARDEN} \
		STATIC=true \
		MIMALLOC=true \
		setup \
		tests \
		install

static-container:
	"${CONTAINER_RUNNER}" build -f "Containerfile" -t ${CONTAINER_NAME} "."

unit-tests-executable: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" unit_tests

unit-tests: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" run-unit-tests

update-regression-tests: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" update-regression-tests

test tests: executables unit-tests regression-tests

.PHONY: clean clean-coverage coverage-xml coverage docs executable \
	executables install regression-tests setup static-container static test \
	tests unit-tests unit-tests-executable update-regression-tests
