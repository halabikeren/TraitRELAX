# TraitRELAX - a tool to associate phenotypic traits with altered selective patterns

TraitRELAX is an open-source software for the joint analysis of binary traits and coding sequence data that allows testing for association of the trait with changes in selection intensity at the codon level across a phylogeny. TraitRELAX is implemented in the C++ library [Bio++](https://github.com/BioPP) (see also: [Bio++ documentation](http://biopp.univ-montp2.fr/)).

## Publication

Halabi K, Levy Karin E, Guéguen L, and Mayrose I. TraitRELAX - A codon model for associating phenotypic traits with altered selective patterns of sequence evolution. Submitted 2020. A preprint is available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.03.04.974584v2)

## Input

The input to TraitRELAX is **a single control file**, which among others, specifies the location of the following files: 
1. A phylogentic tree with branch lengths.
2. A codon multiple sequence alignment (MSA) of the sequence data of the extant species.
3. The character states of the extant species coded as either '0' or '1'.

The TraitRELAX control file specifies parameters as detailed in the [bppSuite manual](http://biopp.univ-montp2.fr/manual/pdf/bppsuite/v0.7.0/bppsuite.pdf). See the provided **TraitRELAX_template.bpp** as an example for such control file.

## Output

TraitRELAX writes the maximum-likelihood scores for the null and alternative models as well as their inferred model parameters to an output file specified in the control file.

## Running the program

Once installed (see next section), the program can be run through the shell:
```
path/to/TraitRELAX/traitrelax param=<path_to_control_file>
```

## Building from source...

### ...by using an installation script

An installation script is available at **install.sh**. To install, modify lines 4,5 to include the desired installation directories of Bio++ and TraitRELAX.

### ...by shell commands
#### Creating source directories
```
bpp_dir=$HOME/local/bpp/
traitrelax_dir=$HOME/local/traitrelax/
mkdir -p $bpp_dir/sources/
mkdir -p $traitrelax_dir
```

#### Downloading the 5 Bio++ libraries
```
cd $bpp_dir/sources
git clone https://github.com/BioPP/bpp-core.git
git clone https://github.com/BioPP/bpp-seq.git
git clone https://github.com/BioPP/bpp-phyl.git
git clone https://github.com/BioPP/bpp-popgen.git
cd $bpp_dir
git clone https://github.com/BioPP/bppsuite.git
```

#### Compiling and installing Bio++
```
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
```

#### Compiling and installing TraitRELAX
```
cd $traitrelax_dir
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .
make -j
make install
```
