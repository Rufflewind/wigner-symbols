#!/bin/sh

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

fcompile dist/tmp/rrf-4.0/src/rrf_module.F90
fcompile dist/tmp/rrf-4.0/src/wigner.F90
flinkshared dist/tmp/lib rrf 4 1 \
    dist/tmp/rrf-4.0/src/rrf_module.F90 \
    dist/tmp/rrf-4.0/src/wigner.F90
