#!/bin/sh
set -euo # "strict" mode
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

PREFIX=$1
TARGET=$2
FORMAT=$3

mkdir -p "${TARGET}"
for exe in adapterremoval3.exe adapterremoval3; do
    if test -e "${PREFIX}/bin/${exe}"; then
        strip "${PREFIX}/bin/${exe}"
        mv -v "${PREFIX}/bin/${exe}" "${TARGET}/"
        break
    fi
done

mv -v "${PREFIX}/share/adapterremoval3/examples" "${TARGET}/examples"
mv -v "${PREFIX}/share/adapterremoval3/docs/html" "${TARGET}/docs"
mv -v "${PREFIX}/share/man/man1/adapterremoval3.1" "${TARGET}/docs/"
cp -v LICENSE "${TARGET}/LICENSE.txt"

if test "${FORMAT}" = "tar.gz"; then
    tar cvzf "${TARGET}.tar.gz" "${TARGET}"
elif test "${FORMAT}" = "zip"; then
    zip -r9 "${TARGET}.zip" "${TARGET}"
else
    echo "Unsupported format: '${FORMAT}'"
    exit 1
fi
