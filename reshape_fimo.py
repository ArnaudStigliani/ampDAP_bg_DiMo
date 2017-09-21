#!/usr/bin/env python2.7

###	This script turns the files we generate after the use of FIMO into the same file shape 
###	than the output of the DNAshapedTFBS algorithm

import sys

def print_predictions(predictions, output):
    """ Print the predictions in the output file. """
    import pandas as pd
    pd_predictions = pd.DataFrame(predictions)
    pd.set_option('display.max_rows', len(pd_predictions))
    with open(output, 'w') as stream:
        stream.write('{0}\n'.format(pd_predictions.to_string(
            index=False, columns=['peak_id', 'start', 'end', 'strand',
                                  'sequence', 'score'])))

def reshape_fimo(fimo_file, fimo_pred):
	with open(fimo_file, 'r') as fimo:
	
		predictions = {'peak_id': [], 'start': [], 'end': [], 'strand': [], 'sequence': [], 'score': []}
		for line in fimo:
			fields = line.rstrip().split()
			predictions['peak_id'].append(fields[2])
			predictions['start'].append(fields[3])
			predictions['end'].append(fields[4])
			predictions['strand'].append(fields[5])
			predictions['sequence'].append(fields[8])
			predictions['score'].append(fields[6])
	
		print_predictions(predictions, fimo_pred)
	
if __name__ == "__main__":
	
	fimo_file = sys.argv[1]
	fimo_pred = sys.argv[2]
	reshape_fimo(fimo_file, fimo_pred)
	
