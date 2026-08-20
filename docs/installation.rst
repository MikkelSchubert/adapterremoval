##############
 Installation
##############

**********************
 Precompiled binaries
**********************

Binaries are provided for 64-bit x86 Linux, ARM-based MacOSX, and Windows at https://github.com/MikkelSchubert/adapterremoval/releases/

************************
 Installing from source
************************

This section describes how to build and install the latest version of AdapterRemoval from source.

AdapterRemoval uses ``make`` to simplify invoking the underlying meson_ build system. Power-users and package managers may invoke ``meson`` directly; use ``make -Bn`` to view the recommended build commands.

Prerequisites
=============

Building and installing AdapterRemoval requires the following:

- ``make``
- A C++17 capable compiler, such as GCC v9+ or Clang v8+. If using GCC, then v11+ is required for AVX512 support
- meson_ v0.56+
- ninja_ v1.10+ (may work with older versions)
- python_ 3.7+
- libdeflate_ v1.7+ (may work with older versions)
- isa-l_ v2.30+
- Sphinx_ is required to build the man page (on by default) and the full documentation (off by default)

Of these dependencies, only libdeflate and isa-l are required to run AdapterRemoval once it has been installed.

- **Debian** and related distributions:

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

    wget -O adapterremoval-3.0.1.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v3.0.1.tar.gz
    tar xvzf adapterremoval-3.0.1.tar.gz
    cd adapterremoval-3.0.1
    make setup
    make

It is furthermore recommended to run the test-suite to verify that AdapterRemoval works correctly on your system:

.. code-block:: console

    make tests

The resulting ``adapterremoval3`` executable is located in the ``build/src`` subdirectory, and can be run as-is. It is also possible to perform a system-wide installation of the AdapterRemoval executable, man-page, and examples using the following command:

.. code-block:: console

    sudo make install

Build options
-------------

The following options may be passed to ``make`` when building AdapterRemoval. For example, to build AdapterRemoval without the man page, thereby removing the dependency on Sphinx:

.. code-block:: console

    make setup MANPAGE=disabled
    make

**Core options:**

- ``MANPAGE=enabled`` / ``MANPAGE=disabled`` / ``MANPAGE=auto`` to enable or disable building of the man page. Requires Sphinx. If set to ``auto``, the man page will only be built if Sphinx is available. Defaults to ``enabled``.
- ``DOCS=enabled`` / ``DOCS=disabled`` / ``DOCS=auto`` to enable or disable building of the full documentation available at https://adapterremoval.readthedocs.org. Requires Sphinx. If set to ``auto``, the documentation page will only be built if Sphinx is available. Defaults to ``disabled``.
- ``STATIC=true`` or ``STATIC=false`` to build statically or dynamically linked binary. Static libraries must be available to build a static binary. Defaults to ``false``.
- ``HARDEN=true`` / ``HARDEN=false`` to enable / disable additional hardening flags during compilation. In some cases it may be necessary to disable these, if the build-system also sets hardening flags. Defaults to ``true``.
- ``PREFIX=/path/to/prefix`` is passed to the meson using the ``--prefix`` option. Is not set by default.

**Development and experimental options:**

- ``COVERAGE=true`` / ``COVERAGE=false`` to include coverage instrumentation in the build. Defaults to ``false``.
- ``DEBUG=true`` / ``DEBUG=false`` to toggle runtime POSIX asserts. Defaults to ``false``, unless ``COVERAGE`` is set to ``true``.
- ``LTO=true`` / ``LTO=false`` to toggle link-time optimization (experimental). Causes irreproducible builds and possibly other issues. Defaults to ``false``.
- ``LTO_MODE=thin`` / ``LTO_MODE=default`` to select the link-time optimization mode. ``thin`` is recommended, but older systems may need to use ``make LTO_MODE=default``. Defaults to ``thin``.
- ``UV=enabled`` / ``UV=disabled`` / ``UV=auto`` to select whether or not uv_ should be used to install optional requirements for running regression tests (fastjsonschema_). Defaults to ``auto``, which uses ``uv >= 0.5.17`` if it is available.
- ``SANITIZE=true`` / ``SANITIZE=false`` to enable address and undefined behavior sanitation at runtime. This has a significant performance overhead. Defaults to ``false``.

********************************
 Building a static Linux binary
********************************

A podman_/docker_ ``Containerfile`` is provided, which is used to generate the pre-compiled binaries mentioned above. To build this, either podman_ or docker_ are required.

To build the container and the static binary, run

.. code-block:: console

    make static-container static

The resulting executable and other files are saved to ``build/static/install``.

.. _docker: https://www.docker.com/

.. _fastjsonschema: https://github.com/horejsek/python-fastjsonschema

.. _homebrew: https://brew.sh

.. _isa-l: https://github.com/intel/isa-l/

.. _libdeflate: https://github.com/ebiggers/libdeflate/

.. _meson: https://mesonbuild.com/

.. _msys2: https://www.msys2.org/

.. _ninja: https://ninja-build.org/

.. _podman: https://podman.io/

.. _python: https://www.python.org/

.. _sphinx: https://www.sphinx-doc.org/

.. _uv: https://docs.astral.sh/uv/
