.. highlight:: Bash

##############
 Installation
##############

********************
 Precompiled binary
********************

Binaries are provided for 64-bit x86 Linux and Windows under https://github.com/MikkelSchubert/adapterremoval/releases/

*************************
 Installing from sources
*************************

This section describes how to build and install the latest version of AdapterRemoval from sources.

Prerequisites
=============

Building and installing AdapterRemoval requires basic build tools including a C++17 capable compiler, meson_ v1.2+, ninja_, python_ 3.8+, libdeflate_ and isa-l_ v2.30+. Sphinx_ is additionally required to build the documentation.

- **Debian**:

  .. code::

     sudo apt-get install build-essential meson ninja-build libdeflate-dev libisal-dev python3 python3-sphinx pkgconf

- **OSX**, requires Homebrew_ to install the dependencies:

  .. code::

     brew install llvm meson ninja isa-l libdeflate sphinx-doc

- **Windows**, requires MSYS2 with an UCRT64 environment:

  .. code::

    pacman -S make mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-isa-l mingw-w64-ucrt-x86_64-libdeflate mingw-w64-ucrt-x86_64-meson mingw-w64-ucrt-x86_64-python mingw-w64-ucrt-x86_64-python-sphinx

Running AdapterRemoval requires only libdeflate and isa-l to be installed.

Building AdapterRemoval
=======================

To compile AdapterRemoval, first download and unpack the newest release from GitHub, and then run ``make``:

.. code::

   wget -O adapterremoval-3.0.0-alpha3.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v3.0.0-alpha3.tar.gz
   tar xvzf adapterremoval-3.0.0-alpha3.tar.gz
   cd adapterremoval-3.0.0-alpha3
   make

The resulting ``adapterremoval3`` executable is located in the 'build/src' subdirectory, and can be run as-is. It is also possible to perform a system-wide installation of the AdapterRemoval executable, man-page, and examples using the following command:

.. code::

   sudo make install

********************************
 Building a static Linux binary
********************************

A podman_/docker_ ``Containerfile`` is provided, which is used to generate the pre-compiled binaries mentioned above. To build this, either podman_ or docker_ are required.

To build the container and the static binary, run

.. code::

   make static-container static

The resulting executable and extra files are saved to ``build/static/install``.

.. _docker: https://www.docker.com/

.. _homebrew: https://brew.sh

.. _isa-l: https://github.com/intel/isa-l/

.. _libdeflate: https://github.com/ebiggers/libdeflate/

.. _meson: https://mesonbuild.com/

.. _ninja: https://ninja-build.org/

.. _podman: https://podman.io/

.. _python: https://www.python.org/

.. _sphinx: https://www.sphinx-doc.org/

