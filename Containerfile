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

ENTRYPOINT [ "/usr/bin/make", "-C", "/root/src", "setup", "test", "regression", "install", "STATIC=true", "BUILDDIR=/root/out/static-container", "DESTDIR=/root/out/static-container/install" ]
