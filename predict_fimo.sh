#!/usr/bin/env bash

meme_motif=$1
bed_input=$2
fasta_input=$3
FIMO_output=$4


threshold=0.5

fimo --text --bgfile --motif-- --verbosity 1 -thresh $threshold $meme_motif \
$fasta_input | sed "1d" | sort -k3,3 -k7,7nr | sort -u -k3,3 > $FIMO_output

FIMO_nb=`wc -l $FIMO_output | cut -f1 -d" " `
bed_nb=`wc -l $bed_input | cut -f1 -d" " `

if [ $FIMO_nb -ne $bed_nb ]; then
	threshold=1
	fimo --text --bgfile --motif-- --verbosity 1 -thresh $threshold $meme_motif \
	$fasta_input | sed "1d" | sort -k3,3 -k7,7nr | sort -u -k3,3 > $FIMO_output
fi
