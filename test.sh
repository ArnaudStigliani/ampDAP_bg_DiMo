#!/usr/bin/env bash


>>../results/stderr.txt
>>../results/stdout.txt

path_scripts=/storage/mathelierarea/processed/arnaudst/ampDAP_bg_DiMO/scripts
data=/storage/mathelierarea/processed/arnaudst/data
results_shared=/storage/mathelierarea/processed/arnaudst/ampDAP_bg_DiMO/results
path_BiasAway=/lsc/BiasAway

genome=$data/Arabidopsis_thaliana.TAIR10.dna.genome.fa
methyl_genome=$data/AraTha_TAIR10_methyl_genome.fa
methyl_alphabet=$data/alphabet
chr_lengths=$data/AraTha_chr_length.genome

#ls -d $data/dap_data_v4/peaks/*/* | grep "colamp" | sed -e 's/_colamp.*//g' | xargs -I{} --max-proc=8 bash -c \

#ls -d $data/dap_data_v4/peaks/C2H2_tnt/At5g66730_colamp_a | grep "colamp" | sed -e 's/_colamp.*//g' | xargs -I{} --max-proc=5 bash -c \
ls -d $data/dap_data_v4/peaks/[N-Z]*/* | grep "colamp" | sed -e 's/_colamp.*//g' | xargs -I{} --max-proc=8 bash -c \
'echo {}

echo "lool"

'
