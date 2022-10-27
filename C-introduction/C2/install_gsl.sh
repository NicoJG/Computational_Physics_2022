#!/bin/bash
# GSL install

die() {
    echo $@
    exit 1
}

install() {
    tmpdir="$(mktemp -d)"
    pushd "$tmpdir"

    wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz || die "Failed downloading archive"
    tar -xvf gsl-2.6.tar.gz || die "Failed extracting archive"
    pushd gsl-2.6/
    ./configure --prefix="$HOME/.local" || die "Failed configure"
    make || die "Failed compilation"
    make install || die "Failed install"

    popd
    popd

    rm -rf "$tmpdir"
}

install

echo "export C_INCLUDE_PATH=\"$HOME/.local/include\"" >> "$HOME/.bashrc"
echo "export LIBRARY_PATH=\"$HOME/.local/lib\"" >> "$HOME/.bashrc"
echo "export LD_LIBRARY_PATH=\"$HOME/.local/lib\"" >> "$HOME/.bashrc"
