##############
 Installation
##############

********************
 Precompiled binary
********************

Binaries are provided for 64-bit x86 Linux, ARM-based MacOSX, and Windows at https://github.com/MikkelSchubert/adapterremoval/releases/

.. warning::

    The static Linux binary uses Musl_ for maximum compatibility, at the cost of a significant performance loss. It is therefore recommended to build AdapterRemoval yourself, or use your distribution's AdapterRemoval package (if available).

************************
 Installing from source
************************

This section describes how to build and install the latest version of AdapterRemoval from source.

Prerequisites
=============

Building and installing AdapterRemoval requires basic build tools including a C++17 capable compiler, meson_ v1.2+, ninja_, python_ 3.9+, libdeflate_ and isa-l_ v2.30+. Sphinx_ is additionally required to build the documentation. Of these dependencies, only libdeflate and isa-l are required to run AdapterRemoval, once it has been installed.

- **Debian**:

.. code-block:: console

    sudo apt-get install build-essential meson ninja-build libdeflate-dev libisal-dev python3 python3-sphinx pkgconf

- **OSX**, requires Homebrew_ to install the dependencies:

.. code-block:: console

    brew install llvm meson ninja isa-l libdeflate sphinx-doc pkgconf

- **Windows**, requires MSYS2_ with an UCRT64 environment:

.. code-block:: console

    pacman -S make mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-isa-l mingw-w64-ucrt-x86_64-libdeflate mingw-w64-ucrt-x86_64-meson mingw-w64-ucrt-x86_64-python mingw-w64-ucrt-x86_64-python-sphinx

Building AdapterRemoval
=======================

To compile AdapterRemoval, first download and unpack the newest release from GitHub, and then run ``make`` in the resulting directory:

.. code-block:: console

    wget -O adapterremoval-3.0.0-alpha3.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v3.0.0-alpha3.tar.gz
    tar xvzf adapterremoval-3.0.0-alpha3.tar.gz
    cd adapterremoval-3.0.0-alpha3
    make

The resulting ``adapterremoval3`` executable is located in the ``build/src`` subdirectory, and can be run as-is. It is also possible to perform a system-wide installation of the AdapterRemoval executable, man-page, and examples using the following command:

.. code-block:: console

    sudo make install

********************************
 Building a static Linux binary
********************************

A podman_/docker_ ``Containerfile`` is provided, which is used to generate the pre-compiled binaries mentioned above. To build this, either podman_ or docker_ are required.

To build the container and the static binary, run

.. code-block:: console

    make static-container static

The resulting executable and other files are saved to ``build/static/install``.

.. _docker: https://www.docker.com/

.. _homebrew: https://brew.sh

.. _isa-l: https://github.com/intel/isa-l/

.. _libdeflate: https://github.com/ebiggers/libdeflate/

.. _meson: https://mesonbuild.com/

.. _msys2: https://www.msys2.org/

.. _musl: https://wiki.musl-libc.org/

.. _ninja: https://ninja-build.org/

.. _podman: https://podman.io/

.. _python: https://www.python.org/

.. _sphinx: https://www.sphinx-doc.org/
