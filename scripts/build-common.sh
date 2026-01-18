#!/bin/sh
set -euo # "strict" mode
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

echo RUNNING SETUP
# WORKAROUND for empty `$@` failing with `set -u` on OSX
if test $# -gt 0; then
    make setup "PREFIX=${PWD}/build/install" DOCS=true "${@}"
else
    make setup "PREFIX=${PWD}/build/install" DOCS=true
fi

echo RUNNING COMPILE
make executables

echo RUNNING TESTS
make tests

echo RUNNING INSTALL
make install

echo RUNNING EXAMPLES
make -C examples EXE="${PWD}/build/src/adapterremoval3"
