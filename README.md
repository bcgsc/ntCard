ntCard 
=
ntCard is a streaming algorithm for cardinality estimation in genomics datasets. As input it takes file(s) in fasta, fastq, sam, or bam formats and computes the total number of distinct k-mers, *F<sub>0</sub>*, and also the *k*-mer coverage frequency histogram, *f<sub>i</sub>*, *i>=1*.  


## Install ntCard on macOS

Install [Homebrew](https://brew.sh/), and run the command

	brew install brewsci/bio/ntcard

## Install ntCard on Linux

Install [Linuxbrew](http://linuxbrew.sh/), and run the command

	brew install brewsci/bio/ntcard

Compiling ntCard from GitHub
===========================

When installing ntCard from GitHub source the following tools are
required:

* [Autoconf](http://www.gnu.org/software/autoconf)
* [Automake](http://www.gnu.org/software/automake)

To generate the configure script and make files:

	./autogen.sh
 
Compiling ntCard from source
===========================
To compile and install ntCard in /usr/local:

```
$ ./configure
$ make 
$ sudo make install 
```

To install ntCard in a specified directory:

```
$ ./configure --prefix=/opt/ntCard
$ make 
$ make install 
```

ntCard uses OpenMP for parallelization, which requires a modern compiler such as GCC 4.2 or greater. If you have an older compiler, it is best to upgrade your compiler if possible. If you have multiple versions of GCC installed, you can specify a different compiler:

```
$ ./configure CC=gcc-xx CXX=g++-xx 
```

For the best performance of ntCard, pass `-O3` flag:  

```
$ ./configure CFLAGS='-g -O3' CXXFLAGS='-g -O3' 
```


To run ntCard, its executables should be found in your PATH. If you installed ntCard in /opt/ntCard, add /opt/ntCard/bin to your PATH:

```
$ PATH=/opt/ntCard/bin:$PATH
```

Run ntCard
==========
```
ntcard [OPTIONS] ... FILE(S) ...
```
Parameters:
  * `-t`,  `--threads=N`: use N parallel threads [1] (N>=2 should be used when input files are >=2)
  * `-k`,  `--kmer=N`: the length of *k*-mer
  * `-c`,  `--cov=N`: the maximum coverage of *k*-mer in output `[1000]`
  * `-p`,  `--pref=STRING`: the prefix for output file names 
  * `-o`,  `--output=STRING`: the name for single output file name (can be used only for single compact output file)
  * `FILE(S)`: input file or set of files seperated by space, in fasta, fastq, sam, and bam formats. The files can also be in compressed (`.gz`, `.bz2`, `.xz`) formats . A list of files containing file names in each row can be passed with `@` prefix.
  
For example to run ntcard on a test file `reads.fastq` with `k=50` and output the histogram in a file with prefix `freq`:
```
$ ntcard -k50 -p freq reads.fastq 
```
To run ntcard on a test file `reads.fastq` with multiple k's `k=32,64,96,128` and output the histograms in files with prefix `freq` use:
```
$ ntcard -k32,64,96,128 -p freq reads.fastq 
```
As another example, to run ntcard on `5` input files file_1.fq.gz, file_2.fa, file_3.sam, file_4.bam, file_5.fq with `k=64` and 6 threads and maximum frequency of `c=100` on output file with prefix `freq`:
```
$ ntcard -k64 -c100 -t6 -p freq file_1.fq.gz file_2.fa file_3.sam file_4.bam file_5.fq
```

If we have a list of input files `lib.in`, to run ntCard with `k=144` and `12` threads and output file with prefix `freq`:
```
$ ntcard -k144 -t12 -p freq @lib.in 
```
Publications
============

## [ntCard](http://bioinformatics.oxfordjournals.org/content/early/2017/01/04/bioinformatics.btw832)

Hamid Mohamadi, Hamza Khan, and Inanc Birol.
**ntCard: a streaming algorithm for cardinality estimation in genomics data**.
*Bioinformatics* (2017) 33 (9): 1324-1330.
[10.1093/bioinformatics/btw832 ](http://dx.doi.org/10.1093/bioinformatics/btw832)
