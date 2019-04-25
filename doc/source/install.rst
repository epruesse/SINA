.. highlight:: shell

Installing SINA
===============

You can install SINA

1. `using Bioconda`_ (recommended)
2. from `pre-compiled tarballs`_ (alternate)
3. or build SINA `from source`_ (for developers)


.. _`using Bioconda`:

Install using Bioconda
----------------------

SINA is available as a Conda_ package in the Bioconda_ channel. Check
the `package info`_ page for more information.

To install, follow these steps:

1. Install Miniconda_ (skip if you've got conda already)

   Download the Miniconda_ installer (links for MacOS_ and Linux_),
   execute it and follow the instructions it shows in the shell::

      # if you are on MacOS
      wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
      # if you are on Linux
      wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

      # then
      sh Miniconda3-lastest-*-x86_64.sh

2. Add the Conda-Forge_ and Bioconda_ channels::

      conda config --add channels defaults
      conda config --add channels bioconda
      conda config --add channels conda-forge

3. Install SINA into a conda environment::

       conda create -n sina sina

4. Activate and use environment::

       conda activate sina
       sina --help


SINA should also work fine if installed with conda into the ``base``
environment. If you encounter problems, try the separate environment.

.. _Conda: https://conda.io
.. _Bioconda: https://bioconda.github.io
.. _Conda-Forge: https://conda-forge.org
.. _Miniconda: https://conda.io/miniconda.html
.. _MacOS: https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
.. _Linux: https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
.. _`package info`: https://bioconda.github.io/recipes/sina/README.html


.. _`pre-compiled tarballs`:

Install from pre-compiled tarballs
----------------------------------

Tar archives containing pre-compiled binaries and requisite libraries
are available on the `SINA releases`_ page at Github. Head on over
there and download the Linux or Macos one. Inside the folder created
by unpacking the archive, you should find a `sina` executable::

  tar xf ~/Downloads/sina-1.6.0-rc.1-linux.tar.gz
  ~/Downloads/sina-1.6.0-rc.1-linux/sina --help

To install SINA system wide, place the contents of the archive in
`/opt` and create symlinks into `/usr/local/bin`::

  wget https://github.com/epruesse/SINA/releases/download/v1.6.0-rc.1/sina-1.5.0-linux.tar.gz
  sudo tar xf sina-1.6.0-rc.1-linux.tar.gz -C /opt
  rm sina-1.6.0-rc.1-linux.tar.gz
  sudo ln -s /opt/sina-1.6.0-rc.1-linux /opt/sina
  sudo ln -s /opt/sina/bin/sina /usr/local/bin/sina

.. _`SINA releases`: https://github.com/epruesse/SINA/releases


.. _`from source`:

Build from source code
----------------------

Building SINA from source can be challenging because SINA depends on
the ARB development libraries. Pre-compiled versions of these are
currently only available from Bioconda, so you may have to start by
building ARB from source.


Prerequisites
~~~~~~~~~~~~~

When building from a source tar ball from the `SINA releases`_ page,
you will need to install:

- ARB_ >= 6.0.0
- Boost_ >= 1.62
  - Boost Thread
  - Boost Program Options
  - Boost IO Streams
  - Boost Filesystem
  - Boost Serialization (SINA < 1.5)
  - Boost System
  - Boost Unit Test Framework (when building / running tests)
- TBB_ >= 2017
- zlib

When building from the raw git source code directly, you will additionally need:

- Autoconf
- Automake
- Libtool
- pkg-config

To build the documentation, you need:

- Sphinx >= 1.8


The easiest way to get all requirements is using Conda::

  if test $(uname) == Linux; then
    dist_extra="gcc patchelf coreutils"
  else
    dist_extra="llvm"
  fi
  conda create -n sina_build automake autoconf libtool pkg-config boost arb-bio-devel \
                  git tbb tbb-devel glib libiconv bc sed sphinx $dist_extra


Building
~~~~~~~~

1. If you are building from git, start by checking out the source code
   and generating the `configure` script::

     git clone https://github.com/epruesse/SINA.git sina
     cd sina
     autoreconf --force --install

2. Run the configure script, pointing it at all required libraries as
   necessary and choosing features and build types::

     ./configure --prefix=install_location \
                 --with-arbhome=path_to_arbhome \
		 --with-boost=path_to_boost_install \
		 --with-boost-libdir=path_to_boost_libs

   If you used conda to install your dependencies, this line should work::

     conda activate sina_build
     mkdir build
     cd build
     ../configure --prefix `pwd`/install \
                  --disable-docs \
		  --with-tbb=$CONDA_PREFIX \
		  --with-boost=$CONDA_PREFIX \
		  --with-boost-libdir=$CONDA_PREFIX/lib \
		  --with-arbhome=$CONDA_PREFIX/lib/arb \
		  LDFLAGS="$LDFLAGS -Wl,-rpath,$CONDA_PREFIX/lib"

   Essential options to **configure**:

   .. program:: configure

   .. option:: --prefix=PATH (/usr/local)

      Set the folder under which ``./bin/sina``, ``./lib/libsina.so`` (or
      ``.dylib``), etc. will be installed.

   .. option:: --with-tbb=PATH

      Set the location of libraries and headers for the Intel
      Threading Building Blocks library.

   .. option:: --with-boost=PATH

      Set the location of the boost header files (without the
      ``include/`` part).

   .. option:: --with-boost-libdir=PATH

      Set the location of the boost lib folder. Often, this is the
      value you used for :option:`--with-boost` with ``/lib``
      appended.

   .. option:: --with-arbhome=PATH

      Set the location of the ARB build directory. Not needed if you
      have ``$ARBHOME`` set to point to the place you built ARB. When
      using the Bioconda package `arb-bio-devel`, use
      `$CONDA_PREFOX/lib/arb` where ``$CONDA_PREFIX`` is the root of
      the environment you installed ARB into.

   .. option:: --with-buildinfo=TEXT

      Set an additional string to be added to the version to identify your build.

   .. option:: --enable-code-coverage

      Add compiler flags to collect code coverage statistics.

   .. option:: --enable-debug

      Enable debug options (sets ``-DDEBUG -O0 -ggdb3`` instead of ``-DNDEBUG -O2 -g``).

   .. option:: --enable-asan

      Enable address sanitizer (sets `-fsanitize=address`).

   .. option:: --enable-fat-tar

      Alters the build so that ``make bindist-gzip`` constructs a fully
      contained tar archive of the build.

   .. option:: --enable-profiling

      Add compiler flags collecting profiling statistics (``-pg``).

   .. option:: --disable-docs

      Do not build the documentation.

   If you installed the dependencies in system wide, standard FHS
   locations, the **configure** script should detect the locations
   correctly. Otherwise you may have to use the ``--with-something``
   options to point it at the right places. If things go wrong, the
   full error messages will be in ``config.log``.

3. Build SINA (replace ``<number of cpus`` with however many cores you've got)::

     make -j<number of cpus>
     make install

   To build binary archives (see also :option:`--enable-fat-tar`), use::

     make bindist-gzip2

   To run unit tests, call::

     make check

   To run only part of the tests, call::

     make check-filtered P=pattern

   where ``pattern`` matches the name(s) of the test you wish to (re)run.

   To run unit tests collecting code coverage, call::

     make check-code-coverage

   To see the full command line of compiler and linker instead of the
   abbreviated display, append ``V=1`` to the ``make`` commandline.



.. _ARB: http://www.arb-home.de
.. _Boost: https://www.boost.org
.. _TBB: https://www.threadingbuildingblocks.org
