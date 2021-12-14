#!/bin/bash

##############################
# crossAlign.sh checks sequencing run samples for their proper
# contents. It aligns the first 10000 reads of each sample to a 
# collection of all possible templates in the run and reports back
# a breakdown of how many reads aligned to which template. It expects
# two arguments: (1) A directory containing the FASTQ files of the run;
# and (2) a folder containing parameters sheets for the DMS experiments
# inside the run. (I will work on expanding this for general purposes)
# 
# Contact: Jochen Weile <jochenweile@gmail.com>
#
# Dependencies: 
#  * binaries: bowtie2, samtools, R
#  * R-packages: yogitools, yogiseq, tileseqMave
#
# usage: ./probe_alignments.sh <FASTQDIR> <PARAMDIR>
#
# output: 
#  * alnProbe.txt : containing a textual representation
#  * visProbe.pdf : A visualization of the results

#helper function to print usage information
usage () {
  cat << EOF

crossAlign.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

Cross-align a sequencing run against possible library to check 
for contamination.
Usage: barseq.sh <FASTQDIR> <FASTA>

<FASTQDIR> : The input directory containing the fastq.gz files
<FASTA>    : A FASTA file containing all possible references

EOF
 exit $1
}


#FASTQ input folder
FQDIR=$1
#libraryFASTA
FASTA=$2

if [[ ! -d $FQDIR ]]; then
  echo "ERROR: First argument must be a directory!" >&2
  usage 1
fi
if [[ ! -r $FASTA || ! $FASTA =~ (.fa|.fasta)$ ]]; then
  echo "ERROR: Second argument must be a FASTA file!">&2
  usage 1
fi
if ! ls ${FQDIR}/*R1*.fastq.gz; then
  echo "ERROR: No FASTQ files found in $FQDIR">&2
  usage 1
fi

#parameter sheet directory
PARAMDIR=$(mktemp -d)

#build bowtie library
bowtie2-build $FASTA $PARAMDIR/probeLibrary

#align 10000 reads from each file to bowtie library and count hits
function runAlignments() {
  #iterate over FASTQ files in input folder
  # echo "Folder: $FQDIR"
  for FQ in $(ls ${FQDIR}/*R1*.fastq.gz); do
    #extract sample ID
    SAMPLE=$(basename $FQ|sed -r "s/_.*//")
    #print sample ID
    echo "Sample $SAMPLE"
    #align the first 40000 lines (=10000 reads), extract reference column
    # then count how often each reference is found and print the result
    bowtie2 --local --very-sensitive-local -x $PARAMDIR/probeLibrary \
      -U <(zcat "$FQ"|head -n 40000) 2>/dev/null\
    |samtools view|awk '{print $3}'|sort|uniq -c 
  done
}
runAlignments >alnProbe.txt

Rscript -e '
options(stringsAsFactors=FALSE)
library(yogiseq)
library(yogitools)

#read input files
infile <- commandArgs(TRUE)[[1]]
fastafile <- commandArgs(TRUE)[[2]]

#read template sequences and extract names
seqs <- yogiseq::readFASTA(fastafile)
templates <- sapply(seqs,function(s)s$getID())
#add "*" for "no match"
templates <- c(templates,"*")

#parse alignment breakdown file
lines <- readLines(infile)
hlines <- grep("^Sample",lines)
starts <- hlines+1
ends <- c(hlines[-1]-1,length(lines))
samples <- yogitools::extract.groups(lines[hlines],"^Sample (.+)$")
#add leading zeros to sample names to allow more intuitive sorting order
samples <- sapply(samples, function(s) {
  if (is.na(as.integer(s))){
    s
  } else {
    sprintf("%03d",as.integer(s))
  }
})

#generate a matrix to store the results
mat <- matrix(
  0, nrow=length(samples),
  ncol=length(templates), 
  dimnames=list(samples,templates)
)
#parse data lines and write them to the matrix
for (i in 1:length(samples)) {
  data <- as.data.frame(yogitools::extract.groups(
    lines[starts[[i]]:ends[[i]]], "^ *(\\d+) +(.+)$"
  ))
  data[,1] <- as.integer(data[,1])/10000
  mat[samples[[i]],data[,2]] <- data[,1]
}

#re-order matrix 
sorder <- order(samples)
mat <- mat[sorder,]
samples <- samples[sorder]

#define color scheme for plot
cmap <- yogitools::colmap(c(0,1),c("white","steelblue3"))
cmapred <- yogitools::colmap(c(0,1),c("white","firebrick3"))
colmat <- apply(mat,1:2,cmap)
colmat[,"*"] <- cmapred(mat[,"*"])

#calculate plotting coordinates
xy <- which(mat > 0,arr.ind=TRUE)
colvals <- colmat[mat>0]

#draw plot
pdf("visProbe.pdf",7,40)
op <- par(las=2,mar=c(10,4,4,1),yaxs="i")
plot(
  NA,type="n",
  xlim=c(0,length(templates)+1),
  ylim=c(c(0,length(samples)+1)),
  axes=FALSE,
  xlab="",ylab="sample"
)
axis(1,1:length(templates),templates)
axis(2,1:length(samples),samples)
rect(xy[,2]-.5,xy[,1]-.5,xy[,2]+.5,xy[,1]+.5,col=colvals,border=NA)
abline(
  h=1:length(samples)-.5,
  v=1:length(templates)-.5,
  col="gray",lty="dotted"
)
par(op)
dev.off()
' alnProbe.txt $FASTA

echo "Done!"

