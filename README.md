### Assignment
This pipeline calculates the coverage of an Illumina MiSeq data set against a reference sequence, 

and tests for correlation between coverage and sequence content.

HIV genome reference: http://www.ncbi.nlm.nih.gov/nuccore/K03455.1

MiSeq data set: http://www.ncbi.nlm.nih.gov/sra/?term=SRR961514

### Commands for running

```
python coverage.py -1 forward.fastq -2 reverse.fastq -x ref.fasta -q 0 -o outdir
```

### Arguments

```
-1  forward fastq 
-2  reverse fastq 
-x  Reference genome fasta
-q  Mean base quality (range 0-41)
```
### Outputs
```
All output will be written into the outdir folder. The default name is "out".

Inside the folder, you will find:

1. coverage.tsv: A tab-delimited file with columns “position” (in the reference sequence)
and “count”. Count is defined as the number of mapped reads that contain the given
reference position in their mapping.

2. coverage.pdf: A plot showing the coverage, with reference position on the x-axis and
count on the y-axis.

3. correlation.html: A short report including tables, plots, or descriptions of statistically 
significant correlation was observed between the nucleotide content of the
reference and the coverage.
```

### References
```
Benjamini, Yuval, and Terence P. Speed. "Summarizing and correcting the GC content bias in 
high-throughput sequencing." Nucleic acids research 40.10 (2012): e72-e72.
```
