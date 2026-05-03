#!/bin/sh
set -euo # "strict" mode
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# WORKAROUND: Custom libdeflate library, since static linking against the
#             corresponding package fails with missing symbols
mkdir -p dependencies/libdeflate

cd dependencies/libdeflate
wget https://github.com/ebiggers/libdeflate/releases/download/v1.25/libdeflate-1.25.tar.gz
tar xvzf libdeflate-1.25.tar.gz

cd libdeflate-1.25

if command -v ccache; then
    cmake -B build \
        -DCMAKE_INSTALL_PREFIX=/ucrt64 \
        -DCMAKE_C_COMPILER_LAUNCHER=ccache \
        -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
else
    cmake -B build \
        -DCMAKE_INSTALL_PREFIX=/ucrt64
fi

cmake --build build
cmake --install build
