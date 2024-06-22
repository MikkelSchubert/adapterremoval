#
# This container is designed to simplify creation of statically linked
# AdapterRemoval executables:
#
#   $ git clone https://github.com/MikkelSchubert/adapterremoval.git
#   $ cd adapterremoval/
#   $ podman build -t ar3builder .
#   $ podman run --rm -it -v${PWD}/:/root/adapterremoval/ ar3builder
#
# This installs AdapterRemoval in build/install/
#
FROM alpine:3.16.2

RUN apk add build-base nasm python3 py3-sphinx cmake

WORKDIR /root


# isa-l - used for compresion (at low gzip levels) and decompression
ARG ISAL=2.30.0
ARG ISAL_HASH=bcf592c04fdfa19e723d2adf53d3e0f4efd5b956bb618fed54a1108d76a6eb56
ADD "https://github.com/intel/isa-l/archive/refs/tags/v${ISAL}.tar.gz" "/root/isa-l-${ISAL}.tar.gz"

RUN echo "${ISAL_HASH}  isa-l-${ISAL}.tar.gz" | sha256sum -c -
RUN tar xzf "isa-l-${ISAL}.tar.gz"
RUN make -C "isa-l-${ISAL}" -j4 -f Makefile.unx install


# libdeflate - used for compression at default gzip levels
ARG DEFLATE=1.14
ARG DEFLATE_HASH=89e7df898c37c3427b0f39aadcf733731321a278771d20fc553f92da8d4808ac
ADD "https://github.com/ebiggers/libdeflate/archive/refs/tags/v${DEFLATE}.tar.gz" "/root/libdeflate-${DEFLATE}.tar.gz"

RUN echo "${DEFLATE_HASH}  libdeflate-${DEFLATE}.tar.gz" | sha256sum -c -
RUN tar xzf "libdeflate-${DEFLATE}.tar.gz"
RUN make -C "libdeflate-${DEFLATE}" -j4 install


# mimalloc - replacement for (slow) musl malloc
ARG MIMALLOC=2.0.7
ARG MIMALLOC_HASH=f23aac6c73594e417af50cb38f1efed88ef1dc14a490f0eff07c7f7b079810a4
ADD "https://github.com/microsoft/mimalloc/archive/refs/tags/v${MIMALLOC}.tar.gz" "/root/mimalloc-${MIMALLOC}.tar.gz"

RUN echo "${MIMALLOC_HASH}  mimalloc-${MIMALLOC}.tar.gz" | sha256sum -c -
RUN tar xzf "mimalloc-${MIMALLOC}.tar.gz"
RUN cmake -S "mimalloc-${MIMALLOC}" -B "mimalloc-${MIMALLOC}/out/release"
RUN make -C "mimalloc-${MIMALLOC}/out/release" -j4 install

ENV LDLIBS="-lmimalloc"
ENV LDFLAGS="-L/usr/local/lib/mimalloc-2.0/"
ENTRYPOINT [ "make", "-C", "/root/adapterremoval/", "CXX=g++", "COLOR=no", "STATIC=yes", "DESTDIR=build/install" ]

CMD ["-j4", "regression", "examples", "test", "install"]
