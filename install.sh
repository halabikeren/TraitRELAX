#! /bin/sh

#### replace with your desired installation directories
bpp_dir=$1
traitrelax_dir=$bpp_dir

#### leave as is

## install Bio++
mkdir -p $bpp_dir/sources/
cd $bpp_dir/sources
git clone https://github.com/BioPP/bpp-core.git
git clone https://github.com/BioPP/bpp-seq.git
git clone -b kerenDevel https://github.com/halabikeren/bpp-phyl.git
cd $bpp_dir/sources/bpp-core
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make
make install
cd ../../bpp-seq
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make
make install
cd ../../bpp-phyl
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make
make install

## install traitrelax
mkdir -p $traitrelax_dir
cd $traitrelax_dir
git clone https://github.com/halabikeren/TraitRELAX.git
cd TraitRELAX
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .
make
make install
