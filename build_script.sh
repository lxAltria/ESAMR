#!/bin/bash

source_dir=`pwd`
external_dir=${source_dir}/external
mkdir -p external
cd ${external_dir}
# build SZ (to use ZSTD compressor)
git clone https://github.com/szcompressor/SZ.git
cd SZ
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_INSTALL_PREFIX=${external_dir}/SZ/install ..
make -j 8
make install

# build SZ3 (to use quantizer and huffman encoder)
cd ${external_dir}
git clone https://github.com/szcompressor/SZ3.git
cp ${source_dir}/SZ3_CMakeLists.txt SZ3/CMakeLists.txt
cd SZ3
mkdir -p build
cd build
cmake ..
make -j 8

# build MGARDx
cd ${external_dir}
git clone https://github.com/lxAltria/MGARDx.git
cd MGARDx
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_INSTALL_PREFIX=${external_dir}/MGARDx/install ..
make -j 8
make install

# build MDR
cd ${source_dir}
mkdir -p build
cd build
cmake ..
make -j 8