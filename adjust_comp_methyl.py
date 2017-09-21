#!/usr/bin/env python2.7

###	This script adjusts the complementary methylated fasta sequences: when using bedtools getfasta with methylated sequences on minus strand, it gives k as complement of m instead of 1 and it lets 1 as complement of itself so we have to change it

import sys
from itertools import izip_longest

def adjust_compl_methyl(methyl_fasta_in, bed_hits, methyl_fasta_out):
	
	strands = {}
	sequences = {}
	
	with open(bed_hits, 'r') as bed:
		for line in bed:
			fields = line.rstrip().split()
			strands[fields[3]] = fields[5]			
			
	with open(methyl_fasta_in, 'r') as fasta_in:
		for line1, line2 in izip_longest(fasta_in, fasta_in, fillvalue=''):
			name = line1.rstrip()
			seq = line2.rstrip()
			if name[0] == ">":
				sequences[name[1:]] = seq
	
	for key in sequences:
		if key in strands:
			if strands[key] == "-":
				adjusted_seq = ""
				for i in range(0,len(sequences[key])):
					if sequences[key][i] == "k":
						adjusted_seq += "1"
						continue
					elif sequences[key][i] == "1":
						adjusted_seq += "m"
						continue
					else:
						adjusted_seq += sequences[key][i]
						continue
				sequences[key] = adjusted_seq
		else:
			print "problem with key {}".format(key)
		
	with open(methyl_fasta_out, 'w') as fasta_out:
		for key in sequences:
			fasta_out.write(">{0}\n{1}\n".format(key, sequences[key]))


if __name__ == "__main__":
	
	methyl_fasta_in = sys.argv[1]
	bed_hits = sys.argv[2]
	methyl_fasta_out = sys.argv[3]
	adjust_compl_methyl(methyl_fasta_in, bed_hits, methyl_fasta_out)
	
