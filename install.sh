#! /bin/sh

#### replace with your desired installation directories
bpp_dir=BIOPP_INSTALLATION_DIRECTORY
traitrelax_dir=TRAITRELAX_INSTALLATION_DIRECTORY

#### leave as is

## install Bio++
mkdir -p $bpp_dir/sources/
cd $bpp_dir
git clone https://github.com/BioPP/bppsuite.git
cd sources
git clone https://github.com/BioPP/bpp-core.git
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make
make install
git clone https://github.com/BioPP/bpp-seq.git
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make
make install
git clone -b kerenDevel https://github.com/halabikeren/bpp-phyl.git
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make
make install
git clone https://github.com/BioPP/bpp-popgen.git
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make
make install
cd ../../../bppsuite/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE ..
make
make install

## install traitrelax
mkdir -p $traitrelax_dir
cd $traitrelax_dir
git clone https://github.com/halabikeren/traitrelax.git
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .
make
make install

