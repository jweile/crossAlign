#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

#dependencies:
library(tileseqMave)
library(yogiseq)
library(argparser)

#process CLI arguments
p <- arg_parser(
  "extract tile sequences from json parameter sheets into a FASTA file",
  name="parajson2fasta.R"
)
p <- add_argument(p, "indir", help="input directory containing json files")
p <- add_argument(p, "--outfile", help="The desired output file name.", default="probeLibrary.fasta")
args <- parse_args(p)

#find json files in the input folder
paramfiles <- list.files(args$indir,pattern="json$",full.names=TRUE)

#function to scrape tile sequences from a parameter sheet json file
scrapeTileSeqs <- function(paramfile) {
  params <- tileseqMave::parseParameters(paramfile)
  geneName <- params$template$geneName
  with(params,setNames(mapply(
    substr,
    x=template$seq,
    start=tiles[,"Start NC in Template"],
    stop=tiles[,"End NC in Template"]
  ),sprintf("%s-Tile%02d",geneName,tiles[,"Tile Number"])))
}

#scrape all input json files for their sequences
sequences <- do.call(c,lapply(paramfiles,scrapeTileSeqs))

#write to output FASTA file
con <- file(args$outfile,open="w")
yogiseq::writeFASTA(con,sequences)
close(con)

cat("Done!\n")
