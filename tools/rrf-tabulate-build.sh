#!/bin/sh
set -eu

fcompile() {
    src=$1
    printf "compiling ... %s\n" "$src"
    f95 -fPIC -O3 -mtune=native -DGFORTRAN -o "$src.o" -c "$src"
}

flinkshared() {
    dir=$1
    name=$2
    major=$3
    rest=$4
    shift 4
    printf "linking ... %s\n" "$dir/lib$name.so"
    mkdir -p "$dir"
    f95 -fPIC -shared "-Wl,-soname,lib$name.so.$major" \
        -o "$dir/lib$name.so.$major.$rest" "$@"
    ln -fs "lib$name.so.$major.$rest" "$dir/lib$name.so.$major"
    ln -fs "lib$name.so.$major" "$dir/lib$name.so"
}

getrrf() (
    dir=$1
    mkdir -p "$dir"
    cd "$dir"
    curl -fLOSs http://www-stone.ch.cam.ac.uk/pub/rrf-4.0.tgz
    cat >rrf-4.0.tgz.sha512sum <<EOF
99feedf949deb9c5576a6d10628e912db4e997053e73914b78f76b04d14783c57026a86eaa91bb50b7240b522ecd59958a2358568a436646f3890aa38b4278b1  rrf-4.0.tgz
EOF
    sha512sum -c rrf-4.0.tgz.sha512sum
    tar xzf rrf-4.0.tgz
)

[ -d dist/tmp/rrf-4.0 ] ||
    getrrf dist/tmp

[ -f dist/tmp/lib/librrf.so ] || {
    fcompile dist/tmp/rrf-4.0/src/rrf_module.F90
    fcompile dist/tmp/rrf-4.0/src/wigner.F90
    flinkshared dist/tmp/lib rrf 4 1 \
        dist/tmp/rrf-4.0/src/rrf_module.F90 \
        dist/tmp/rrf-4.0/src/wigner.F90
}

mkdir -p dist/tmp/bin
cabal exec -- \
    ghc -Wall -O -o dist/tmp/bin/rrf-tabulate tools/rrf-tabulate.hs \
        -Ldist/tmp/lib -lrrf
LD_LIBRARY_PATH=dist/tmp/lib dist/tmp/bin/rrf-tabulate
md5sum dist/rrf*.txt
