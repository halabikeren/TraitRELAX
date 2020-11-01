# TraitRELAX - a tool to associate phenotypic traits with altered selective patterns

[ ![Codeship Status for halabikeren/TraitRELAX](https://app.codeship.com/projects/1108cd20-6cd5-0138-bfdb-1e3b1ab831af/status?branch=master)](https://app.codeship.com/projects/394727)

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

TraitRELAX writes the maximum-likelihood scores for the null and alternative models as well as their inferred model parameters to STDOUT. You can save the results by redirecting STDOUT into a file (see Examples/README.txt). 
Additionaly it can save results to output files specified in the control file (see the provided **TraitRELAX_template.bpp**).

## Running the program

Once installed (see next section), the program can be run through the shell:
```
path/to/TraitRELAX/traitrelax param=<path_to_control_file>
```

## Installing via docker

Rather than installing the program from scratch, you can pull a docker image with the pre-compiled program. To do so, first install [Docker](https://docs.docker.com/get-docker/) on your machine and then run on the command line:
```
docker pull halabikeren/traitrelax:version1.0.0
```
To run the program, first create on your machine a directory with the input for traitrelax, including the input data and a control file (see the Examples folder for more details). Then, run the following on the command line:
```
docker run -v <path_to_input_directory>:/traitrelax/exec/ -it traitrelax param=<name_of_control_file>
```

## Building from source...

### ...by using an installation script

An installation script is available at **install.sh**. To install, run on the command line:
```
sh install.sh <path_to_prgram>
```
Once the installation is complete, the progeam will be available in ```<path_to_prgram>/TraitRELAX/TraitRELAX/TraitRELAX```

### ...by shell commands

The compilation may take a little while (especially that of `bpp-phyl`; ~20-30 min using 16-cores) so perhaps make some tea
#### Creating source directories
```
bpp_dir=BIOPP_INSTALLATION_DIRECTORY
mkdir -p $bpp_dir/sources/
```

#### Downloading the 4 Bio++ libraries
```
cd $bpp_dir/sources
git clone https://github.com/BioPP/bpp-core.git
git clone https://github.com/BioPP/bpp-seq.git
git clone -b kerenDevel https://github.com/halabikeren/bpp-phyl.git
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
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE OMP_NUM_THREADS=<required_number_of_cores> ..
make -j
make install
```

#### Compiling and installing TraitRELAX
clone
```
traitrelax_dir=TRAITRELAX_INSTALLATION_DIRECTORY # the directory to which you clone TraitRELAX
mkdir -p $traitrelax_dir
cd $traitrelax_dir/
git clone https://github.com/halabikeren/TraitRELAX.git
```
and compile
```
cd $traitrelax_dir/TraitRELAX/
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE .
make -j
make install
```

#### Running TraitRELAX with Parallelization

To run TraitRELAX with multiple threads, you should state the number of required cores in the cmake command of the complication of the Bio++ bpp-phyl library:
```
OMP_NUM_THREADS=<required_number_of_cores>
```

#### Distinction from the RELAX method

TraitRELAX varies from The RELAX basic method in several aspects:

1. TraitRELAX does not require a prior partition of the branches into two categories.
2. TraitRELAX uses a fixed effect likelihood approach for the selective categories, and as such constrains the selective category of each site to remain consistent across branches. That being said, the selection value of a site can change upon relaxation or intensification based on the inferred evolutionary history of the examined trait.
3. TraitRELAX uses a different optimization approach that limits the search space of the selection intenisty parameter to (0,10] in aim of avoiding extreme estimates.

Due to these differences, estimation of the parameters shared between the two methods may be inconsistent with one another.
