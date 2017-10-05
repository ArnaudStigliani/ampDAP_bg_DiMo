#!/usr/bin/env bash

path_scripts=../scripts
results_data=../results/DAP_DNAshape
results_shared=../results


motif_input=$1
meme_matrix=$2
dimo_matrix=$3
methylation=$4


#Turn motif into matrix for DiMO input
description_line=`grep -n "letter-probability matrix" $motif_input | awk -F [:] '{print $1}'`
line_count=`wc -l $motif_input | awk '{print $1}'`

tail -n $(($line_count-$description_line)) $motif_input > $meme_matrix

#Convert matrix from MEME motif to matrix for DiMO input
python $path_scripts/convert_matrix.py $meme_matrix $dimo_matrix $methylation

#Add name of motif
sed -i '1i > motif_input' $dimo_matrix
