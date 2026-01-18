#!/bin/sh
set -euo # "strict" mode

pacman -S --noconfirm \
    make \
    mingw-w64-ucrt-x86_64-gcc \
    mingw-w64-ucrt-x86_64-isa-l  \
    mingw-w64-ucrt-x86_64-libdeflate \
    mingw-w64-ucrt-x86_64-meson \
    mingw-w64-ucrt-x86_64-mimalloc \
    mingw-w64-ucrt-x86_64-python \
    mingw-w64-ucrt-x86_64-python-jsonschema \
    mingw-w64-ucrt-x86_64-python-sphinx \
    mingw-w64-ucrt-x86_64-uv
