#!/bin/sh
set -euo # "strict" mode

STATIC=0
while test $# -gt 0; do
	case "${1}" in
		--static)
			STATIC=1
			echo "Enabling static build"
			;;
		*)
			error "ERROR: Unknown option '${1}'"
			exit 1
			;;
	esac

	shift 1
done

sudo apt-get install -y \
	build-essential \
	libdeflate-dev \
	meson \
	ninja-build \
	pkgconf \
	python3 \
	python3-fastjsonschema \
	python3-sphinx

# libisal-dev does not include static library
if test ${STATIC} -eq 1; then
	sudo apt-get install -y \
		nasm \
        autoconf \
        libtool

	export ISAL_V="2.32.0"
	wget "https://github.com/intel/isa-l/archive/refs/tags/v${ISAL_V}.tar.gz" -O isa-l-${ISAL_V}.tar.gz
	echo "7a194ff80d0f7e20615c497654e8a51b0184d0c79e2e265c7f555f52a26a05a4  isa-l-${ISAL_V}.tar.gz" | sha256sum -c -
	tar xvzf isa-l-${ISAL_V}.tar.gz
	cd isa-l-${ISAL_V}

    export DEB_BUILD_MAINT_OPTIONS="hardening=+all"
    eval $(dpkg-buildflags --export=sh)
    ./autogen.sh
    ./configure

	make checks
	make tests
	sudo make install V=1
else
	sudo apt-get install -y \
		libisal-dev
fi

