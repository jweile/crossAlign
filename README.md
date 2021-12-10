# Cross-align

Cross-align a sequencing run against possible library to check 
for contamination.

## crossAlign.sh usage:
```bash
crossAlign.sh <FASTQDIR> <FASTA>
```

, where `<FASTQDIR>` is a directory containing fastq.gz files, and 
`<FASTA>` is a FASTA file containing all expected reference sequences. 
A FASTA file with the PhiX genome can be found in `phix.fasta`.

Another tool, `parajson2fasta.R` is provided to automatically extract 
tile sequences from tileseq parameter sheets into a FASTA file.

## parajson2fasta.R usage:
```bash
parajson2fasta.R [-o <OUTFILE>] <JSONDIR>
```