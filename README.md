[![Release](https://img.shields.io/github/release/bcgsc/ntCard.svg)](https://github.com/bcgsc/ntCard/releases)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/ntCard/total?logo=github)](https://github.com/bcgsc/ntCard/archive/master.zip)
[![Conda](https://img.shields.io/conda/dn/bioconda/ntcard?label=Conda)](https://anaconda.org/bioconda/ntcard)
[![Issues](https://img.shields.io/github/issues/bcgsc/ntCard.svg)](https://github.com/bcgsc/ntCard/issues)

![Logo](./ntcard-logo.png)

ntCard is a streaming algorithm for cardinality estimation in genomics datasets. Input any number of file(s) in FASTA, FASTQ, SAM, or BAM formats, and ntCard will output the total number of distinct k-mers ($F_0$) and the *k*-mer coverage frequency histogram ($f_i; i \geq 1$).

# Compiling from source

Make sure you have the following dependencies installed:

- [Meson](https://mesonbuild.com/) >= 1.0.0
- [btllib](https://github.com/bcgsc/btllib) - if installing from source, run btllib's `./compile` script and add the following environment variables:
```
export CPPFLAGS="-isystem /path/to/btllib/install/include $CPPFLAGS"
export LDFLAGS="-L/path/to/btllib/install/lib -lbtllib $LDFLAGS"
```

Clone the repo using `git clone --recurse-submodules https://github.com/bcgsc/ntCard.git`.

Run the following commands in the project's directory to install ntCard:

```
meson setup --buildtype release --prefix=$(pwd)/install build
cd build
ninja install
```

Verify that the installation was successful by executing the `ntcard` binary in `install/bin`.

You can change the installation prefix by modifying the `--prefix` argument of the `meson setup` command. If the `--prefix` argument is removed from the command, ntCard will be installed to the system's default path for binaries.

# Running ntCard

```
Usage: ntcard [options] files 

Estimates k-mer coverage histogram in input files

Positional arguments:
files                   [required]

Optional arguments:
-h --help               shows help message and exits [default: false]
-v --version            prints version information and exits [default: false]
-k --kmer-length        Length of the k-mers, ignored if using a spaced seed (see -s)
-s --spaced-seed        Spaced seed pattern with 1's as cares and 0's as don't cares
-t --num-threads        Number of threads [default: 1]
-c --max-coverage       Maximum coverage to estimate [required]
-l --left-bits          Number of bits to take from the left for sampling [default: 7]
-r --right-bits         Number of bits to take from the right as k-mer representations [default: 27]
--long-mode             Optimize file reader for long sequences (>5kb) [default: false]
```

# Publications

## [ntCard](http://bioinformatics.oxfordjournals.org/content/early/2017/01/04/bioinformatics.btw832)

Hamid Mohamadi, Hamza Khan, and Inanc Birol.
**ntCard: a streaming algorithm for cardinality estimation in genomics data**.
*Bioinformatics* (2017) 33 (9): 1324-1330.
[10.1093/bioinformatics/btw832 ](http://dx.doi.org/10.1093/bioinformatics/btw832)
