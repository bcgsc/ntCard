ntCard 
=
ntCard is a streaming algorithm for cardinality estimation in genomics datasets. As iput it takes file(s) is fastq, fastq, sam, or bam formats and computes the total number of distinct k-mers, *F<sub>0</sub>*, and also the *k*-mer coverage frequency histogram, *f<sub>i</sub>*, *i>=1*.  

# Build the binary for ntcard
Run:
```
$ make
```
# Run ntcard
```
ntcard [OPTIONS] ... [FILE]
```
Parameters:
  * `-k`,  `--kmer=SIZE`: the length of *k*-mer `[64]`
  * `-t`,  `--threads=N`: use N parallel threads `[1]`
  * `-c`,  `--cov=N`: the maximum coverage of *k*-mer in output `[64]`
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

