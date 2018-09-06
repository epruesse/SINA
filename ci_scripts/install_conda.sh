#!/bin/bash
#
# Installs Bio-Conda and all dependencies
# - adds miniconda PATH to bashrc (for CircleCI)
# - manages a tarball if file `LOCAL` present (for ./cicleci build)
# - creates conda_state.txt for caching with circleci "{{checksum}}"
# - updates packages with every run
# - needs BASH_ENV to point to the bashrc
# - needs MINICONDA to point to the miniconda install path


CONDA_PACKAGES="autoconf automake libtool toolchain pkg-config boost arb-bio-devel"
CONDA_PACKAGES="$CONDA_PACKAGES tbb"

case "$(uname)" in
    Linux)
	CONDA_OSNAME=Linux
	CONDA_PACKAGES="$CONDA_PACKAGES gcc patchelf"
	;;
    Darwin)
	CONDA_OSNAME=MacOSX
	CONDA_PACKAGES="$CONDA_PACKAGES llvm"
	;;
esac

# Make sure we have conda in the PATH always
if test -n "$BASH_ENV"; then
    echo "Prepending $MINICONDA/bin to PATH in $BASH_ENV"
    echo export PATH="$MINICONDA/bin:$PATH" >> $BASH_ENV
fi
export PATH="$MINICONDA/bin:$PATH" >> $BASH_ENV

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



