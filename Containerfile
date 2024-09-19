#
# This container is designed to simplify creation of statically linked
# AdapterRemoval executables. See `make static-container` and `make static`.
#
FROM alpine:3.20.1

RUN apk add \
    build-base \
    isa-l-dev \
    isa-l-static \
    libdeflate-dev \
    libdeflate-static \
    meson \
    mimalloc2-dev \
    ninja \
    pkgconf \
    py3-sphinx \
    python3

WORKDIR /root

ENV LDFLAGS="-lmimalloc"
ENV OUT_DIR="/root/out/static-container"
ENV SRC_DIR="/root/src"

CMD meson setup --reconfigure "${OUT_DIR}/builddir" "${SRC_DIR}" -Dstatic=true \
    && meson compile -v -C "${OUT_DIR}/builddir" \
    && meson test --print-errorlogs -C "${OUT_DIR}/builddir" \
    && DESTDIR="${OUT_DIR}/install" meson install -C "${OUT_DIR}/builddir"
