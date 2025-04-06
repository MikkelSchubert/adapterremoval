###############################################################################
# This makefile is provided as a simple convenience

## Optional features; set to 'true' to enable or 'false' to disable:
##   $ make DEBUG=true

# Include coverage instrumentation in build
COVERAGE := false

# Debug build; adds warnings, debugging symbols, and extra STL/POSIX asserts
DEBUG := ${COVERAGE}

# Enable address and undefined behavior sanitation
SANITIZE := false

# Enable hardening flags
HARDEN := false

# Generate statically linked binary
# It is recommended to use the included Containerfile to build the static binary
STATIC := false

# Use/require mimalloc. Mainly intended for static builds, as the musl allocator
# comes with a significant performance cost
MIMALLOC := $(STATIC)

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
	rm -rf "${BUILDDIR}"

coverage: ${NINJAFILE}
	ninja -C "${BUILDDIR}" coverage-text
	cat build/meson-logs/coverage.txt

coverage-xml: ${NINJAFILE}
	ninja -C "${BUILDDIR}" coverage-xml

install: ${NINJAFILE}
	meson install -C "${BUILDDIR}"

regression-tests: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" run-regression-tests

setup ${NINJAFILE}:
	meson setup "${BUILDDIR}" --reconfigure \
		-Db_coverage=${COVERAGE} \
		-Ddebug=${DEBUG} \
		-Dharden=${HARDEN} \
		-Dmimalloc=${MIMALLOC} \
		-Dstatic=${STATIC} \
		${MESON_OPTIONS}

static: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" static

static-container: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" static-container

unit-tests-executable: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" unit_tests

unit-tests: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" run-unit-tests

update-regression-tests: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" update-regression-tests

tests: executables unit-tests regression-tests

.PHONY: clean coverage-xml coverage executable executables install \
	regression-tests setup static-container static unit-tests \
	unit-tests-executable update-regression-tests tests
