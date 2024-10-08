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

###############################################################################

ifeq ($(strip ${DEBUG}), true)
BUILD_OPTS := debugoptimized
else
BUILD_OPTS := release
endif

ifeq ($(strip ${SANITIZE}), true)
SANITIZE_OPTS := address,undefined
else
SANITIZE_OPTS := none
endif

###############################################################################

.PHONY := clean executable install regression setup static static-container test
# Meson commands cannot be run in parallel
.NOTPARALLEL:

# Default meson build directory
BUILDDIR := builddir
# Location of ninja build-file; used to detected existing setup
NINJAFILE := ${BUILDDIR}/build.ninja

executable: ${NINJAFILE}
	meson compile -C "${BUILDDIR}"

clean:
	rm -rvf "${BUILDDIR}"

install: ${NINJAFILE}
	meson install -C "${BUILDDIR}"

regression: ${NINJAFILE}
	meson test -C "${BUILDDIR}" --suite regression

setup:
	meson setup "${BUILDDIR}" --reconfigure \
		--buildtype=${BUILD_OPTS} \
		-Db_coverage=${COVERAGE} \
		-Db_sanitize=${SANITIZE_OPTS} \
		-Dharden=${HARDEN} \
		-Dstatic=${STATIC}

static: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" static

static-container: ${NINJAFILE}
	meson compile -C "${BUILDDIR}" static-container

test: setup
	meson test -C "${BUILDDIR}" --suite unit

${NINJAFILE}:
	$(MAKE) setup
