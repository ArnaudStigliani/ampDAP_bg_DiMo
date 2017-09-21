#!/usr/bin/env python2.7

import os
import sys
import numpy as np


def meme2DiMO(matrix_meme, matrix_DiMO, methylation=False):
	"""
	This function convert a MEME formatted count matrix to a matrix suitable for DiMO input. 
	It first checks if the matrix corresponds to the methylation state declared
	"""	
	lines_in = []
		
	with open(matrix_meme, 'r') as matrix_in:
		for line in matrix_in:
			if line.strip():
				lines_in.append(line.rstrip().split())
	
	try:
		if methylation and len(lines_in[0]) == 6:	
			alphabet = ["A", "C", "G", "T", "m", "1"]
		elif not methylation and len(lines_in[0]) == 4:
			alphabet = ["A", "C", "G", "T"]
		else:
			raise Exception
			
	except Exception as e:
		print "The input matrix length does not correspond to the methylation state you provided. It should be 6 with methylation or 4 whithout methylation, it is not the case"
	
	else:
		array_in = np.array(lines_in)
		array_out = np.transpose(array_in)
		
		with open(matrix_DiMO, 'w') as matrix_out:
			for i in range(0, len(array_out)):
				current_line = "\t".join(str(frequency) for frequency in array_out[i])
				matrix_out.write("{0} |\t{1}\n".format(alphabet[i], current_line))


if __name__ == "__main__":
	
	matrix_meme = sys.argv[1]
	matrix_DiMO = sys.argv[2]
	methylation_state = sys.argv[3]
	
	if methylation_state in ["true", "True", "T", "t", "1"]:
		methylation = True
	else:
		methylation = False

	meme2DiMO(matrix_meme, matrix_DiMO, methylation)
