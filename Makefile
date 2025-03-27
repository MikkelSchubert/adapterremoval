###############################################################################
# This makefile is provided as a simple convenience

## Optional features; set to 'true' to enable or 'false' to disable

# Debug build; adds warnings and debugging symbols
DEBUG := false

# Include coverage instrumentation in build
COVERAGE := false

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

ifeq ($(strip ${DEBUG}), true)
BUILD := debugoptimized
else
ifeq ($(strip ${DEBUG}), false)
BUILD := release
else
$(error "DEBUG must be 'true' or 'false', not '${DEBUG}'")
endif
endif

ifeq ($(strip ${SANITIZE}), true)
SANITIZE_OPTS := address,undefined
else
ifeq ($(strip ${SANITIZE}), false)
SANITIZE_OPTS := none
else
$(error "SANITIZE must be 'true' or 'false', not '${SANITIZE}'")
endif
endif

###############################################################################

.PHONY := clean executable install regression setup static static-container test
# Meson commands cannot be run in parallel
.NOTPARALLEL:

# Default meson build directory
BUILDDIR := build
# Location of ninja build-file; used to detected existing setup
NINJAFILE := ${BUILDDIR}/build.ninja

executable: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" adapterremoval3

clean:
	rm -rvf "${BUILDDIR}"

coverage: ${NINJAFILE}
	ninja -C "${BUILDDIR}" coverage-text
	cat build/meson-logs/coverage.txt

install: ${NINJAFILE}
	meson install -C "${BUILDDIR}"

regression: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" regression

setup:
	meson setup "${BUILDDIR}" --reconfigure \
		--buildtype=${BUILD} \
		-Db_coverage=${COVERAGE} \
		-Db_sanitize=${SANITIZE_OPTS} \
		-Dharden=${HARDEN} \
		-Dstatic=${STATIC} \
		-Dmimalloc=${MIMALLOC}

static: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" static

static-container: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" static-container

test_executable: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" unit_tests

test: ${NINJAFILE}
	meson test -C "${BUILDDIR}" --print-errorlogs --suite unit

${NINJAFILE}:
	$(MAKE) setup
