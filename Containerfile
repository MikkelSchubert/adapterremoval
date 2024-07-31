#
# This container is designed to simplify creation of statically linked
# AdapterRemoval executables. See `make static-container` and `make static`.
#
FROM alpine:3.20.1

RUN apk add build-base nasm python3 py3-sphinx cmake

WORKDIR /root


# isa-l - used for compresion (at low gzip levels) and decompression
ARG ISAL=2.31.0
ARG ISAL_HASH=e218b7b2e241cfb8e8b68f54a6e5eed80968cc387c4b1af03708b54e9fb236f1
ADD "https://github.com/intel/isa-l/archive/refs/tags/v${ISAL}.tar.gz" "/root/isa-l-${ISAL}.tar.gz"

RUN echo "${ISAL_HASH}  isa-l-${ISAL}.tar.gz" | sha256sum -c -
RUN tar xzf "isa-l-${ISAL}.tar.gz"
RUN make -C "isa-l-${ISAL}" -j4 -f Makefile.unx install


# libdeflate - used for compression at default gzip levels
ARG DEFLATE=1.20
ARG DEFLATE_HASH=ed1454166ced78913ff3809870a4005b7170a6fd30767dc478a09b96847b9c2a
ADD "https://github.com/ebiggers/libdeflate/archive/refs/tags/v${DEFLATE}.tar.gz" "/root/libdeflate-${DEFLATE}.tar.gz"

RUN echo "${DEFLATE_HASH}  libdeflate-${DEFLATE}.tar.gz" | sha256sum -c -
RUN tar xzf "libdeflate-${DEFLATE}.tar.gz"
RUN cd "libdeflate-${DEFLATE}" && cmake -B build && cmake --build build && cmake --install build


# mimalloc - replacement for (slow) musl malloc
ARG MIMALLOC=2.1.7
ARG MIMALLOC_HASH=0eed39319f139afde8515010ff59baf24de9e47ea316a315398e8027d198202d
ADD "https://github.com/microsoft/mimalloc/archive/refs/tags/v${MIMALLOC}.tar.gz" "/root/mimalloc-${MIMALLOC}.tar.gz"

RUN echo "${MIMALLOC_HASH}  mimalloc-${MIMALLOC}.tar.gz" | sha256sum -c -
RUN tar xzf "mimalloc-${MIMALLOC}.tar.gz"
RUN cmake -S "mimalloc-${MIMALLOC}" -B "mimalloc-${MIMALLOC}/out/release"
RUN make -C "mimalloc-${MIMALLOC}/out/release" -j4 install

ENV LDLIBS="-lmimalloc"
ENV LDFLAGS="-L/usr/local/lib/mimalloc-2.1/"
ENTRYPOINT [ "make", "-C", "/root/adapterremoval/" ]

CMD ["-j4", "regression", "examples", "test", "install"]
