.. highlight:: Bash

Installation
============

Precompiled binary
------------------

A pre-compiled binary is provided for 64-bit Linux systems under
https://github.com/MikkelSchubert/adapterremoval/releases/


Installing from sources
-----------------------

Installing AdapterRemoval from sources requires `libdeflate`_ and `isa-l`_. On Debian based systems, these may be installed as follows::

    sudo apt-get install libdeflate-dev libisal-dev

In addtion, a C++11 compatible compiler and basic build-tools are required. On Debian based systems, these may be installed as follows::

    sudo apt-get install build-essential

To compile AdapterRemoval, first download and unpack the newest release from GitHub, and then run the 'make' command::

    wget -O adapterremoval-3.0.0-alpha1.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v3.0.0-alpha1.tar.gz
    tar xvzf adapterremoval-3.0.0-alpha1.tar.gz
    cd adapterremoval-3.0.0-alpha1
    make

The resulting 'adapterremoval3' executable is located in the 'build/release' subdirectory, and can be run as-is. It is also possible to perform a system-wide installation of the AdapterRemoval executable, man-page, and examples using the following command::

    sudo make install

.. _libdeflate: https://github.com/ebiggers/libdeflate/
.. _isa-l: https://github.com/intel/isa-l/
