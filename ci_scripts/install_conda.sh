#!/bin/bash

CONDA_PACKAGES="autoconf automake libtool toolchain pkg-config boost gcc arb-bio"
CONDA_PACKAGES="$CONDA_PACKAGES pcre libiconv"

# Make sure we have conda in the PATH always
if ! grep $MINICONDA/bin $BASH_ENVx; then
   echo "Prepending $MINICONDA/bin to PATH in $BASH_ENV"
   echo export PATH="$MINICONDA/bin:$PATH" >> $BASH_ENV
   source $BASH_ENV
fi

# Homemade restore_cache
if test -e "LOCAL"; then
    echo "Running locally"
    if test -e conda.tgz; then
	echo -n "Unpacking local conda cache... "
	tar -C / -xzf conda.tgz
	echo "done"
    fi
fi

# Install Miniconda if missing
if test -d $MINICONDA; then
    echo "Found conda install"
else
    case `uname` in
	Linux)  CONDA_OSNAME=Linux;;
	Darwin) CONDA_OSNAME=MacOSX;;
    esac

    CONDA_BASEURL=https://repo.continuum.io/miniconda
    curl $CONDA_BASEURL/Miniconda3-latest-$CONDA_OSNAME-x86_64.sh -o miniconda.sh

    bash miniconda.sh -b -p $MINICONDA
    hash -r

    conda config --system --set always_yes yes --set changeps1 no
    conda update -q conda
    for channel in r defaults conda-forge bioconda; do
	conda config --system --add channels $channel
    done
    conda update -q conda
fi

# Install/update package
# (Hope that "install" suffices)
conda install -q conda $CONDA_PACKAGES
conda update -q conda $CONDA_PACKAGES
conda info

# Homemade save_cache
# Clean out tarballs, update cache if there were any (=> something was installed)
if test -e "LOCAL"; then
    if conda clean --yes --tarballs | grep Removed; then
	echo -n "Packing local conda cache..."
	tar -C / -czf conda.tgz $MINICONDA
	echo "done"
    fi
else
    conda clean --yes --tarballs
    (conda info; conda list) >> conda_state.txt
fi



