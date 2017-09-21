#!/usr/bin/env python2.7

###	This script removes any sequence containing ambiguous bases (not ACGT) from a list of input sequences
### It takes as input the names of two fasta files, the first one containing the input sequences (potentially with ambiguous bases) and the second one for the output sequences (i.e. all the sequences containing only ACGT) 


import sys
from itertools import izip_longest


def remove_ambiguous(fasta_in, fasta_out):
	DNA_alphabet = 'ACGT'

	with open(fasta_in, 'r') as in_file:
		with open(fasta_out, 'w') as out_file:
			for line1, line2 in izip_longest(in_file, in_file, fillvalue=''):
				name = line1.rstrip()
				seq = line2.rstrip()
				if all(base in DNA_alphabet for base in seq):
					out_file.write('{0}\n{1}\n'.format(name, seq))
				else:
					print "Sequence {0} {1} has been removed".format(name, seq)


if __name__ == "__main__":
	
	fasta_in = sys.argv[1]
	fasta_out = sys.argv[2]
	remove_ambiguous(fasta_in, fasta_out)
