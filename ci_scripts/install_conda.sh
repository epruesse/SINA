#!/bin/bash

case `uname` in
    Linux)  CONDA_OSNAME=Linux;;
    Darwin) CONDA_OSNAME=MacOSX;;
esac

CONDA_BASEURL=https://repo.continuum.io/miniconda
wget $CONDA_BASEURL/Miniconda3-latest-$CONDA_OSNAME-x86_64.sh -O miniconda.sh

bash miniconda.sh -b -p $HOME/miniconda

export PATH="$HOME/miniconda/bin:$PATH"
export PREFIX="$HOME/miniconda"
hash -r

conda config --set always_yes yes --set changeps1 no
conda update -q conda
for channel in r defaults conda-forge bioconda; do
    conda config --system --add channels $channel
done

conda create -n build arb-bio boost pkg-config autoconf automake libtool
source activate build

export CFLAGS="$CFLAGS -I$CONDA_PREFIX/include"
export LDFLAGS="$LDFLAGS -L$CONDA_PREFIX/lib -Wl,-rpath -Wl,$CONDA_PREFIX/lib"
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:$CONDA_PREFIX/share/pkgconfig"
