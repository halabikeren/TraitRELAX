# TraitRELAX - Association of phenotypes with selection intensity at the codon level across a phylogeny

TraitRELAX is an open-source software for joint analysis of binary traits and coding sequence data that allows testing for association of the trait with changes in selection intensity at the codon level across a phylogeny. TraitRELAX is implemented in [Bio++](https://github.com/BioPP). For more information, click [here](http://biopp.univ-montp2.fr/).

## Building from Source

#### Installing from bash

Installation script is available at install.sh. To install, modify lines 4,5 to include the desired installation directories of Bio++ and TraitRELAX.

#### Create source directories

`bpp_dir=$HOME/local/bpp/`
`traitrelax_dir=$HOME/local/traitrelax/`
`mkdir -p $bpp_dir/sources/`
`mkdir -p $traitrelax_dir`

#### Download all Bio++ libraries

`cd $bpp_dir/sources`

`git clone https://github.com/BioPP/bpp-core.git`

`git clone https://github.com/BioPP/bpp-seq.git`

`git clone -b kerenDevel https://github.com/halabikeren/bpp-phyl.git`

`git clone https://github.com/BioPP/bpp-popgen.git`

`cd $bpp_dir`

`git clone https://github.com/BioPP/bppsuite.git`


#### Compile and install

`cd $bpp_dir/sources/bpp-core`

`mkdir build`

`cd build`

`cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .. # prepare compilation`

`make # compile`

`make install # move files to the installation directory`


`cd bpp-seq`

`mkdir build`

`cd build`

`cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .. # prepare compilation`

`make # compile`

`make install # move files to the installation directory`


`cd ../../bpp-phyl`

`mkdir build`

`cd build`

`cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .. # prepare compilation`

`make # compile`

`make install # move files to the installation directory`


`cd ../../bpp-popgen`

`mkdir build`

`cd build`

`cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .. # prepare compilation`

`make # compile`

`make install # move files to the installation directory`


`cd ../../../bppsuite`

`mkdir build`

`cd build`

`cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .. # prepare compilation`

`make # compile`

`make install # move files to the installation directory`


#### Run via command line
`$traitrelax_dir/TraitRELAX/traitrelax param=<path_to_input_parameters_file>`  
+ _`<path_to_input_parameters_file>` is the full path to a parameters file build based on [bppSuite manual](http://biopp.univ-montp2.fr/manual/pdf/bppsuite/v0.7.0/bppsuite.pdf). An example file is available above (see TraitRELAX_template.bpp).
or 
