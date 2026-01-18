#!/bin/sh
set -euo # "strict" mode

sudo apt-get install -y \
	build-essential \
	libdeflate-dev \
	libisal-dev \
	meson \
	ninja-build \
	pkgconf \
	python3 \
	python3-jsonschema \
	python3-sphinx
