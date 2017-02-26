ntCard 
=
ntCard is a streaming algorithm for cardinality estimation in genomics datasets. As iput it takes file(s) is fastq, fastq, sam, or bam formats and computes the total number of distinct k-mers, *F<sub>0</sub>*, and also the *k*-mer coverage frequency histogram, *f<sub>i</sub>*, *i>=1*.  


## Install ntCard on Mac OS X

Install [Homebrew](http://brew.sh/), and run the command

	brew install homebrew/science/ntcard


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
ntcard [OPTIONS] ... [FILE]
```
Parameters:
  * `-k`,  `--kmer=SIZE`: the length of *k*-mer `[64]`
  * `-t`,  `--threads=N`: use N parallel threads `[1]`
  * `-c`,  `--cov=N`: the maximum coverage of *k*-mer in output `[64]`
  * `-p`,  `--pref=STRING`: the prefix for output file name `[freq]`
  * `FILE`: input file or set of files seperated by space, in fasta, fastq, sam, and bam formats. The files can also be in compressed (`.gz`, `.bz2`, `.xz`) formats . A list of files containing file names in each row can be passed with `@` prefix.
  
For example to run ntcard on a test file `reads.fastq` with `k=50`:
```
$ ntcard -k50 reads.fastq 
```
To run ntcard on a test file `reads.fastq` with multiple k's `k=32,64,96,128` use:
```
$ ntcard -k32,64,96,128 reads.fastq 
```
As another example, to run ntcard on `5` input files file_1.fq.gz, file_2.fa, file_3.sam, file_4.bam, file_5.fq with `k=64` and 6 threads and maximum output of frequencies `c=100`:
```
$ ntcard -k64 -c100 -t6 file_1.fq.gz file_2.fa file_3.sam file_4.bam file_5.fq
```

If we have a list of input files `lib.in` with input file names in each row and want to run ntCard with `k=144` and 12 threads:
```
$ ntcard -k144 -t12 @lib.in 
```
Publications
============

## [ntCard](http://bioinformatics.oxfordjournals.org/content/early/2017/01/04/bioinformatics.btw832)

Hamid Mohamadi, Hamza Khan, and Inanc Birol.
**ntCard: a streaming algorithm for cardinality estimation in genomics data**.
*Bioinformatics* (2017).
[10.1093/bioinformatics/btw832 ](http://dx.doi.org/10.1093/bioinformatics/btw832)
