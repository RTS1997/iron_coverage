# iron_coverage
Nextflow script to check the coverage of ont reads.

1. Takes a read.fastq.gz file and a assembled genome assembly.fa file.
2. It uses minimap2 to align the reads to the assembled genome. This outputs a SAM file with the alignment.
3. This alignment is then sorted using samtools, and a sorted BAM file is produced.
4. From this the 'samtools depth' command is used to calculate the 'average read depth' and 'standard deviation' of this alignment.

### Run the program: 

```
./nextflow code/2023_projects/iron_coverage/scripts/ont_coverage_calc.nf --gen_file <fasta file location> --reads_file <compressed reads file location> --output_fld <results file output folder> -resume
```