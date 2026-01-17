#!/bin/sh
set -euo # "strict" mode

echo RUNNING SETUP
# WORKAROUND for empty `$@` failing with `set -u` on OSX
if test $# -gt 0; then
    make setup DEBUG=true SANITIZE=true "${@}"
else
    make setup DEBUG=true SANITIZE=true
fi

echo RUNNING COMPILE
make executables

echo RUNNING TESTS
make tests

echo RUNNING INSTALL
make install DESTDIR="${PWD}/install"

echo RUNNING EXAMPLES
make -C examples EXE="${PWD}/build/src/adapterremoval3"
