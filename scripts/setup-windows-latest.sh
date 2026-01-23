#!/bin/sh
set -euo # "strict" mode
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

pacman -S --noconfirm \
    git \
    make \
    mingw-w64-ucrt-x86_64-cmake \
    mingw-w64-ucrt-x86_64-gcc \
    mingw-w64-ucrt-x86_64-isa-l  \
    mingw-w64-ucrt-x86_64-meson \
    mingw-w64-ucrt-x86_64-mimalloc \
    mingw-w64-ucrt-x86_64-python \
    mingw-w64-ucrt-x86_64-python-jsonschema \
    mingw-w64-ucrt-x86_64-python-sphinx \
    mingw-w64-ucrt-x86_64-uv


# WORKAROUND: Custom libdeflate library, since static linking against the
#             corresponding package fails with missing symbols
mkdir -p dependencies/libdeflate

cd dependencies/libdeflate
wget https://github.com/ebiggers/libdeflate/releases/download/v1.25/libdeflate-1.25.tar.gz
tar xvzf libdeflate-1.25.tar.gz

cd libdeflate-1.25
cmake -B build -DCMAKE_INSTALL_PREFIX=/ucrt64
cmake --build build
cmake --install build
