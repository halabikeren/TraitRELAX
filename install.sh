#! /bin/sh

#### replace with your desired installation directories
bpp_dir=BIOPP_INSTALLATION_DIRECTORY
traitrelax_dir=TRAITRELAX_INSTALLATION_DIRECTORY

#### leave as is

## install Bio++
mkdir -p $bpp_dir/sources/
cd $bpp_dir/sources
git clone https://github.com/BioPP/bpp-core.git
git clone https://github.com/BioPP/bpp-seq.git
git clone https://github.com/BioPP/bpp-phyl.git
git clone https://github.com/BioPP/bpp-popgen.git
cd $bpp_dir
git clone https://github.com/BioPP/bppsuite.git
cd $bpp_dir/sources/bpp-core
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make -j
make install
cd ../../bpp-seq
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make -j
make install
cd ../../bpp-phyl
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make -j
make install
cd ../../bpp-popgen
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make -j
make install
cd ../../../bppsuite
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make -j
make install

## install traitrelax
mkdir -p $traitrelax_dir
cd $traitrelax_dir
git clone https://github.com/halabikeren/traitrelax.git
cd TraitRELAX
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .
make -j
make install

