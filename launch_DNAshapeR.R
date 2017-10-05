#!/usr/bin/Rscript

### This R script uses the DNAshapeR module from "T.-P. Chiu , F. Comoglio, T. Zhou, L. Yang, R. Paro, and R. Rohs: DNAshapeR: an R/Bioconductor package for DNA shape prediction and feature encoding. Bioinformatics 32, 1211-1213 (2016)."
### It calculates different feature values regarding the DNA shapes: Helix Twist (HelT), Propeller Twist (ProT), Minor Groove Width (MGW) and Roll
### As input, it takes 3 arguments + a list of shapes. The arguments are the fasta sequences file on which the values will be calculated, the file containing the methylation positions (as a non conventional fasta file with one line containing the usual header (>sequence_header) and another line containing the CpG methylation positions instead of the usual sequence (pos1, pos2...)) and the name of the output file. The shape list can be a combination of the different shapes cited above (HelT, ProT, MGW and Roll)
###The output will then contain lines of coma separeted values for each shape of the list and for each position of each sequence (one line per sequence)


working_dir <- ".."
setwd(working_dir)

suppressMessages(library(DNAshapeR, lib.loc = "~/.R/"))

args <- commandArgs(trailingOnly=TRUE)

fasta_file <- args[1]
methyl_file <- args[2]
out_file <- args[3]
in_shapes=args[-c(1:3)]

#Calculating all the shape features from the input file (w/ or w/o methylation)
pred <- getShape(fasta_file, methylate = TRUE, methylatedPosFile = methyl_file)

shapes <- vector("character", length(in_shapes))

#Adaptation of the names from DNAshapedTFBS to DNAshapeR
for(i in 1:length(in_shapes)){
	shapes[i] <- switch(in_shapes[i], 
		"HelT" = "1-HelT", 
		"ProT" = "1-ProT",
		"MGW" = "1-MGW",
		"Roll" = "1-Roll")
}

#Retrieving the shape features
featureVector <- encodeSeqShape(fasta_file, pred, shapes)

pasted_features <- vector("character", length(featureVector[,1]))
for(i in 1:length(pasted_features)) {
	pasted_features[i] <- paste(featureVector[i,], collapse=',')
}
write.table(pasted_features, file=out_file, quote=FALSE, row.names=FALSE, col.names=FALSE)

