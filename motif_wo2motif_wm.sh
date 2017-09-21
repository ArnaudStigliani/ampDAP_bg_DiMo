#!/usr/bin/env bash

data=/storage/mathelierarea/processed/Victor/data

genome=$data/Arabidopsis_thaliana.TAIR10.dna.genome.fa
methyl_genome=$data/AraTha_TAIR10_methyl_genome.fa
alphabet=$data/alphabet

### Methylation boolean: if methylation is true, it means that you want to get the methylated motif, starting from the unmethylated one, and it works the other way if the "methylation" boolean is false. Default is false.
methylation=false

OPTIND=1

while getopts mf:b:h:s:g:o: flag; do
	case "${flag}"
	in
		m) methylation=true;;
		f) FIMO_out=$OPTARG;;
		b) bed_hits_file=$OPTARG;;
		h) fasta_hits_file=$OPTARG;;
		s) sites_file=$OPTARG;;
		g) background=$OPTARG;;
		o) out_motif=$OPTARG;;
		\?) echo "Wrong usage of the parameters"; exit 1;;
		*) echo "Wrong usage of the parameters"; exit 1;;
	esac
done



cut -f3-6 $FIMO_out | awk -v OFS='\t' '{split($1, a, "[:-]"); $(($2+=a[2]-1)); $(($3+=a[2])); print a[1],$2,$3,$1,0,$4}' > $bed_hits_file

if [ $methylation = true ]; then
	bedtools getfasta -s -name -fi $methyl_genome -bed $bed_hits_file -fo $fasta_hits_file
else
	bedtools getfasta -s -name -fi $genome -bed $bed_hits_file -fo $fasta_hits_file

sed '/>/d' $fasta_hits_file > $sites_file

if [ $methylation = true ]; then
	sites2meme -ext .sites -alph $alphabet -bg $background > $out_motif
else
	sites2meme -ext .sites -bg $background > $out_motif


