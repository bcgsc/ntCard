ntCard 
=
ntCard is a streaming algorithm for cardinality estimation in genomics datasets.

# Build the binary
Run:
```
$ make
```
# Run ntcard
```
ntcard [OPTIONS] ... [FILE]
```
Parameters:
  * `-k`,  `--kmer=SIZE`: the length of k-mer used for runtime test hashing `[64]`
  * `-t`,  `--threads=N`: use N parallel threads `[1]`
  * `-c`,  `--cov=N`: the maximum coverage of kmer in output `[64]`
  * `FILE`: is the input file or set of files sepreted by space in fasta, fastq, sam, bam formats. The files can also be in compressed (.gz, .bz2, .xz) formats . A list of files containing file names in each rows can be passe with @ prefix.
  
For example to run ntcard on a test file `reads.fastq` with `k=50`:
```
$ ntcard -k50 reads.fastq 
```
