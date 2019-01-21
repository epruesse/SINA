#!/bin/bash
#
# Installs Bio-Conda and all dependencies
# - adds miniconda PATH to bashrc (for CircleCI)
# - manages a tarball if file `LOCAL` present (for ./cicleci build)
# - creates conda_state.txt for caching with circleci "{{checksum}}"
# - updates packages with every run
# - needs BASH_ENV to point to the bashrc
# - needs MINICONDA to point to the miniconda install path
set -x

CONDA_PACKAGES="autoconf automake libtool toolchain pkg-config boost arb-bio-devel"
CONDA_PACKAGES="$CONDA_PACKAGES git tbb tbb-devel glib libiconv bc sed sphinx"

CONDA_BASEURL=https://repo.continuum.io/miniconda

# expand '~' in MINICONDA path (alternatives to eval are too long)
eval MINICONDA=$MINICONDA
export MINICONDA

case "$(uname)" in
    Linux)
	CONDA_OSNAME=Linux
	CONDA_PACKAGES="$CONDA_PACKAGES gxx_linux-64 patchelf coreutils"
	;;
    Darwin)
	CONDA_OSNAME=MacOSX
	CONDA_PACKAGES="$CONDA_PACKAGES clangxx_osx-64"
	;;
esac

# Setup PATH
if test -n "$BASH_ENV"; then
    echo "Prepending $MINICONDA/bin to PATH in $BASH_ENV"
    echo export PATH="$MINICONDA/bin:$PATH" >> $BASH_ENV
fi
export PATH="$MINICONDA/bin:$PATH" >> $BASH_ENV

# Install Miniconda if missing
if test -d $MINICONDA; then
    echo "Found conda install"
else
    curl $CONDA_BASEURL/Miniconda3-latest-$CONDA_OSNAME-x86_64.sh -o miniconda.sh
    bash miniconda.sh -b -p $MINICONDA
    hash -r

    conda config --system --set always_yes yes --set changeps1 no
    conda config --system --add channels defaults
    conda config --system --add channels bioconda
    conda config --system --add channels conda-forge
    conda update -q conda
fi

# Install/update package
conda install -q conda $CONDA_PACKAGES
conda update -q conda $CONDA_PACKAGES
conda info
conda clean --yes --all

# Dump status
mkdir -p conda
conda info > conda/info.txt
conda list > conda/root.txt
ls -1 $MINICONDA/pkgs > conda/pkgs.txt

## Fix for as yet incomplete removal of .la files by conda-forge
find $MINICONDA/lib -name \*.la -delete



