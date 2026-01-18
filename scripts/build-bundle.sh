#!/bin/sh
set -euo # "strict" mode
trap 's=$?; echo >&2 "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

PREFIX=$1
TARGET=$2

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


# WORKAROUND: libdeflate cannot be linked statically
if test -e "/ucrt64/bin/libdeflate.dll"; then
    cp "/ucrt64/bin/libdeflate.dll" "${TARGET}"
    cp "/ucrt64/share/licenses/libdeflate/COPYING" "${TARGET}/docs/libdeflate.COPYING"
fi

tar cvzf "${TARGET}.tar.gz" "${TARGET}"
