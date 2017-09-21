#!/usr/bin/env python2.7

### This script creates a methylated genome from the normal genome and methylation positions. It takes as input a fasta file containing the whole normal genome, a tsv file containing the methylation information and the name of the output methylated genome as fasta file.
### The methylation information file should have the same configuration than the tsv file of methylated positions from Arabidopsis thaliana Col-0 leaf given in "R. Schmitz et al., Patterns of population epigenomic diversity, Nature 495, 193–198, 14 March 2013". (doi:10.1038/nature11968, GEO accession number: GSM1085222)
### The alphabet used to represent the covalent changes is the one from "C. Viner et al., Modeling methyl-sensitive transcription factor motifs with an expanded epigenetic alphabet" posted on March 16, 2016 on BioRxiv (doi: https://doi.org/10.1101/043794), except that we use only the methylation changes (m if a cytosine is methylated and 1 for methylation of the complementary cytosine of a guanine)
 

import sys
import time
from operator import itemgetter

positions = []

start_time = time.time()
with open('/storage/scratch/victor/DNAshape/data/GSM1085222_mC_calls_Col_0.tsv', 'r') as methyl_in:
	methyl_in.readline()
	for line in methyl_in:
		features = line.rstrip().split()
		if int(features[6]) == 1:
			positions.append({'chr':int(features[0]), 'pos':int(features[1]), 'strand':features[2]})
			
positions = sorted(positions, key=itemgetter('chr', 'pos'))

runtime = time.time() - start_time
print 'Time retrieving methylation positions: {:.5f}'.format(runtime)
print 'Memory used is {}'.format(sys.getsizeof(positions))

start_time = time.time()
with open('/storage/scratch/victor/DNAshape/data/Arabidopsis_thaliana.TAIR10.dna.genome.fa', 'r') as ori_genome:
	seq = [[],[],[],[],[]]
	for line in ori_genome:
		if line[0]==">":
			ref_chr = int(line.rstrip()[-1])			#To modify if using human genome (doesn't apply to chr X and Y)
		else:
			for char in line.rstrip():
				seq[ref_chr-1].append(char)

runtime = time.time() - start_time
print 'Time retrieving genome sequences: {:.5f}'.format(runtime)
print 'Memory used is approximatively {}'.format(len(seq)*sys.getsizeof(seq[0]))

start_time = time.time()
for el in positions:
	if el['strand']=='+':
		if seq[el['chr']-1][el['pos']] == 'C':
			seq[el['chr']-1][el['pos']] = 'm'
		else:
			print "Problem: base at position {0} on chrom {1} is {2} instead of C".format(el['pos'], el['chr'], seq[el['chr']-1][el['pos']])
	else:
		if seq[el['chr']-1][el['pos']] == 'G':
			seq[el['chr']-1][el['pos']] = '1'
		else:
			print "Problem: base at position {0} on chrom {1} is {2} instead of G".format(el['pos'], el['chr'], seq[el['chr']-1][el['pos']])

runtime = time.time() - start_time
print 'Time replacing methylation: {:.5f}'.format(runtime)
print 'Memory used is {}'.format(len(seq)*sys.getsizeof(seq[0]))

start_time = time.time()
with open('/storage/scratch/victor/DNAshape/data/AraTha_TAIR10_methyl_genome.fa', 'w') as new_genome:
	for item in seq:
		new_genome.write('>chr{}'.format(seq.index(item)+1))
		i = 0
		while i < len(item):
			if i % 60 == 0:
				new_genome.write('\n')
			new_genome.write(item[i])
			i+=1
		new_genome.write('\n')	

runtime = time.time() - start_time
print 'Time writing methylated genome file: {:.5f}'.format(runtime)
