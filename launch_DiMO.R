#!/usr/bin/Rscript

### This R script uses DiMO (Discriminative Motif Optimizer) from "R.Y. Patel, G.D. Stormo, Discriminative motif optimization based on perceptron training, Bioinformatics, 30 (2014), pp. 941-948"
### If methylation is considered, we use a modified version of DiMO where the tool has been adapted to methylation by adding the "m" and "1" coding charachters (which represent methylation on a cytosine and on the paired cytosine of a guanine respectively) to the non ambiguous DNA alphabet "ACGT".
### As input, it takes 5 arguments: a boolean indicating whether methylation is considered, the positive fasta sequences, the negative fasta sequences, the motif to optimize and the name of the output files
### The output will then contain the optimized motif matrix


working_dir <- ".."
setwd(working_dir)

args <- commandArgs(trailingOnly=TRUE)

methylation <- args[1]
foreground_seq <- args[2]
background_seq <- args[3]
matrix_in <- args[4]
output_name <- args[5]

if (eval(methylation)){
  library(DiMO.methyl.AUPRC, lib.loc = "~/.R/")
} else {
  library(DiMO.AUPRC, lib.loc = "~/.R/")
}

DiMO(
  positive.file=foreground_seq,
  negative.file=background_seq,
  pfm.file.name=matrix_in,
  output.flag=output_name,
  epochs=150,
  add.at.five=0,
  add.at.three=0,
  learning.rates=seq(1,0.1,length.out=3)
)
