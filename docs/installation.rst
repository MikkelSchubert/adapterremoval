.. highlight:: Bash

Installation
============


Installation with Conda
-----------------------

If you have `Conda`_ installed on your system::

    conda install -c bioconda adapterremoval


Installing on Debian based systems
----------------------------------

Debian users on Stretch, Buster, or Sid, or using Jessie-backports, as well as Ubuntu users on Zesty or Artful, may install AdapterRemoval using apt::

	apt-get install adapterremoval

For other distributions, or to get the latest version of AdapteRemoval, please see the `Installing from sources`_ section below.


Installing on OSX
-----------------

MacOSX users may install AdapterRemoval using Homebrew::

	brew install homebrew/science/adapterremoval

Please see the Homebrew website for instructions on how to install and use Homebrew:

    https://brew.sh/


Installing from sources
-----------------------

Installing AdapterRemoval from sources requires the presence of libz and bz2 headers. On Debian based systems, these may be installed as follows::

    sudo apt-get install zlib1g-dev libbz2-dev

In addtion, a C++11 compatible compiler and basic build-tools are required. On Debian based systems, these may be installed as follows::

    sudo apt-get install build-essential

To compile AdapterRemoval, first download and unpack the newest release from GitHub, and then run the 'make' command::

    wget -O adapterremoval-2.3.4.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v2.3.4.tar.gz
    tar xvzf adapterremoval-2.3.4.tar.gz
    cd adapterremoval-2.3.4
    make

The resulting 'AdapterRemoval' executable is located in the 'build' subdirectory, and can be run as-is. It is also possible to perform a system-wide installation of the AdapterRemoval executable, man-page, and examples using the following command::

    sudo make install


.. _Conda: https://conda.io/docs/