#!/bin/sh
set -euo # "strict" mode

echo RUNNING SETUP
make setup DEBUG=true SANITIZE=true "${@:-}"

echo RUNNING COMPILE
make executables

echo RUNNING TESTS
make tests

echo RUNNING INSTALL
make install DESTDIR="${PWD}/install"

echo RUNNING EXAMPLES
make -C examples EXE="${PWD}/build/src/adapterremoval3"
