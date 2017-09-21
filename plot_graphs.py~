#!/usr/bin/env python2.7

import os
import sys
from sklearn import metrics as skm
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_methyl_sensitivity(sensitivity_file):
	"""
	Function to parse the methylation sensitivity from the input file. The input file contains 3 columns, tab separated, with the family of the TF in the first column, its name in the second one and the sensitivity (float) in the last one. The file should a first line with the names of the columns. 
	It returns a dictionary with the TF reference as keys (family/name). The values of the dictionary are lists with 13 elements, the first one for the sensitivity of the TF, the other 12 for the scores (AUPRC and AUROC) obtained with the different methods 
	"""
	scores = {}
	
	with open(sensitivity_file, 'r') as sensitivity:
		sensitivity.readline()
		for line in sensitivity:
			features = line.rstrip().split()
			
			TF = "{0}/{1}".format(features[0], features[1])
			
			scores[TF] = [float(features[2]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
			
	return scores


def parse_scores(AUPRC_file, AUROC_file, methyl_sensitivity_file):
	"""
	Function to parse the AUPRC and AUROC scores in the input files. Each scores file has 4 tab separated columns: family of the TF, name of the TF, method + methylation state used for the predictions and the score obtained. See the function "parse_methyl_sensitivity" for more information about the methylation sensitivity file. Only the TFs for which we have information about their methylation sensitivity are considered.
	This function returns a list of tuples with the following structure: (TF_family/TF_name, [methylation sensitivity, 6 AUPRC scores, 6 AUROC scores]). That list is sorted by increasing methylation sensitivity.
	"""
	scores = parse_methyl_sensitivity(methyl_sensitivity_file)
	final_scores = {}
	
	with open(AUPRC_file, 'r') as AUPRC_data:
		AUPRC_data.readline()
		for line in AUPRC_data:
			features = line.rstrip().split()
			TF = '{0}/{1}'.format(features[0], features[1])
			if TF in scores:										#Only the TFs with methylation sensitivity data
				if features[2] == "from_methyl":
					if features[3] == "wm_FIMO_wm_shapes":
						scores[TF][1] = float(features[4])
					elif features[3] == "wo_FIMO_wm_shapes":
						scores[TF][2] = float(features[4])
					elif features[3] == "wm_FIMO_wo_shapes":
						scores[TF][3] = float(features[4])
					elif features[3] == "wo_FIMO_wo_shapes":
						scores[TF][4] = float(features[4])
					elif features[3] == "wm_FIMO":
						scores[TF][5] = float(features[4])
					elif features[3] == "wo_FIMO":
						scores[TF][6] = float(features[4])
					else:
						print "Problem with line {0} in file {1}".format(line, AUPRC_file)
				else:
					if features[3] == "wm_FIMO_wm_shapes":
						scores[TF][7] = float(features[4])
					elif features[3] == "wo_FIMO_wm_shapes":
						scores[TF][8] = float(features[4])
					elif features[3] == "wm_FIMO_wo_shapes":
						scores[TF][9] = float(features[4])
					elif features[3] == "wo_FIMO_wo_shapes":
						scores[TF][10] = float(features[4])
					elif features[3] == "wm_FIMO":
						scores[TF][11] = float(features[4])
					elif features[3] == "wo_FIMO":
						scores[TF][12] = float(features[4])
					else:
						print "Problem with line {0} in file {1}".format(line, AUPRC_file)
					
			
	with open(AUROC_file, 'r') as AUROC_data:
		AUROC_data.readline()
		for line in AUROC_data:
			features = line.rstrip().split()
			TF = '{0}/{1}'.format(features[0], features[1])
			if TF in scores:										#Only the TFs with methylation sensitivity data
				if features[2] == "from_methyl":
					if features[3] == "wm_FIMO_wm_shapes":
						scores[TF][13] = float(features[4])
					elif features[3] == "wo_FIMO_wm_shapes":
						scores[TF][14] = float(features[4])
					elif features[3] == "wm_FIMO_wo_shapes":
						scores[TF][15] = float(features[4])
					elif features[3] == "wo_FIMO_wo_shapes":
						scores[TF][16] = float(features[4])
					elif features[3] == "wm_FIMO":
						scores[TF][17] = float(features[4])
					elif features[3] == "wo_FIMO":
						scores[TF][18] = float(features[4])
					else:
						print "Problem with line {0} in file {1}".format(line, AUROC_file)
				else:
					if features[3] == "wm_FIMO_wm_shapes":
						scores[TF][19] = float(features[4])
					elif features[3] == "wo_FIMO_wm_shapes":
						scores[TF][20] = float(features[4])
					elif features[3] == "wm_FIMO_wo_shapes":
						scores[TF][21] = float(features[4])
					elif features[3] == "wo_FIMO_wo_shapes":
						scores[TF][22] = float(features[4])
					elif features[3] == "wm_FIMO":
						scores[TF][23] = float(features[4])
					elif features[3] == "wo_FIMO":
						scores[TF][24] = float(features[4])
					else:
						print "Problem with line {0} in file {1}".format(line, AUROC_file)
	
	
	#Remove elements for which methylation sensitivity is available but no score has been calculated
	for TF in scores:
		if scores[TF][1] != 0.0:
			final_scores[TF] = scores[TF]
		
	#Sort by methylation sensitivity, then by AUPRC score from methyl motif and FIMO w/ methyl data and DNAshape w/ methyl data)
	scores_sorted = sorted(final_scores.iteritems(), key=lambda x: (x[1][0], x[1][1]))
	
	return scores_sorted


def classify_by_sensitivity(scores):
	"""
	This function takes all the scores and for each type of score, it creates a dictionary with 3 lists to classify by methylation sensitivity: methyl_plus for the TFs that bind preferentially when cytosines are methylated, methyl_minus for the TFs that bind preferentially when cytosines are NOT methylated and non_sensitive when the TF binding is not affected by methylation.
	The function returns these 12 dictionaries.
	"""
	AUPRC_from_methyl_wm_FIMO_wm_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUPRC_from_methyl_wo_FIMO_wm_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUPRC_from_methyl_wm_FIMO_wo_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUPRC_from_methyl_wo_FIMO_wo_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUPRC_from_methyl_wm_FIMO = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUPRC_from_methyl_wo_FIMO = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	
	AUPRC_from_non_methyl_wm_FIMO_wm_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUPRC_from_non_methyl_wo_FIMO_wm_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUPRC_from_non_methyl_wm_FIMO_wo_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUPRC_from_non_methyl_wo_FIMO_wo_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUPRC_from_non_methyl_wm_FIMO = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUPRC_from_non_methyl_wo_FIMO = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	
	
	AUROC_from_methyl_wm_FIMO_wm_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUROC_from_methyl_wo_FIMO_wm_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUROC_from_methyl_wm_FIMO_wo_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUROC_from_methyl_wo_FIMO_wo_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUROC_from_methyl_wm_FIMO = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUROC_from_methyl_wo_FIMO = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	
	AUROC_from_non_methyl_wm_FIMO_wm_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUROC_from_non_methyl_wo_FIMO_wm_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUROC_from_non_methyl_wm_FIMO_wo_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUROC_from_non_methyl_wo_FIMO_wo_shapes = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUROC_from_non_methyl_wm_FIMO = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	AUROC_from_non_methyl_wo_FIMO = {'methyl_plus':[], 'non_sensitive':[], 'methyl_minus':[]}
	
#	labels = []
	
	for TF in scores:
#		labels.append(TF[0])
		if TF[1][0] >= 1.0:
			AUPRC_from_methyl_wm_FIMO_wm_shapes['methyl_plus'].append(TF[1][1])
			AUPRC_from_methyl_wm_FIMO_wm_shapes['non_sensitive'].append(np.nan)
			AUPRC_from_methyl_wm_FIMO_wm_shapes['methyl_minus'].append(np.nan)
			
			AUPRC_from_methyl_wo_FIMO_wm_shapes['methyl_plus'].append(TF[1][2])
			AUPRC_from_methyl_wo_FIMO_wm_shapes['non_sensitive'].append(np.nan)
			AUPRC_from_methyl_wo_FIMO_wm_shapes['methyl_minus'].append(np.nan)			
		
			AUPRC_from_methyl_wm_FIMO_wo_shapes['methyl_plus'].append(TF[1][3])
			AUPRC_from_methyl_wm_FIMO_wo_shapes['non_sensitive'].append(np.nan)
			AUPRC_from_methyl_wm_FIMO_wo_shapes['methyl_minus'].append(np.nan)
			
			AUPRC_from_methyl_wo_FIMO_wo_shapes['methyl_plus'].append(TF[1][4])
			AUPRC_from_methyl_wo_FIMO_wo_shapes['non_sensitive'].append(np.nan)
			AUPRC_from_methyl_wo_FIMO_wo_shapes['methyl_minus'].append(np.nan)
			
			AUPRC_from_methyl_wm_FIMO['methyl_plus'].append(TF[1][5])
			AUPRC_from_methyl_wm_FIMO['non_sensitive'].append(np.nan)
			AUPRC_from_methyl_wm_FIMO['methyl_minus'].append(np.nan)
			
			AUPRC_from_methyl_wo_FIMO['methyl_plus'].append(TF[1][6])
			AUPRC_from_methyl_wo_FIMO['non_sensitive'].append(np.nan)
			AUPRC_from_methyl_wo_FIMO['methyl_minus'].append(np.nan)

			
			AUROC_from_methyl_wm_FIMO_wm_shapes['methyl_plus'].append(TF[1][7])
			AUROC_from_methyl_wm_FIMO_wm_shapes['non_sensitive'].append(np.nan)
			AUROC_from_methyl_wm_FIMO_wm_shapes['methyl_minus'].append(np.nan)
			
			AUROC_from_methyl_wo_FIMO_wm_shapes['methyl_plus'].append(TF[1][8])
			AUROC_from_methyl_wo_FIMO_wm_shapes['non_sensitive'].append(np.nan)
			AUROC_from_methyl_wo_FIMO_wm_shapes['methyl_minus'].append(np.nan)			
		
			AUROC_from_methyl_wm_FIMO_wo_shapes['methyl_plus'].append(TF[1][9])
			AUROC_from_methyl_wm_FIMO_wo_shapes['non_sensitive'].append(np.nan)
			AUROC_from_methyl_wm_FIMO_wo_shapes['methyl_minus'].append(np.nan)
			
			AUROC_from_methyl_wo_FIMO_wo_shapes['methyl_plus'].append(TF[1][10])
			AUROC_from_methyl_wo_FIMO_wo_shapes['non_sensitive'].append(np.nan)
			AUROC_from_methyl_wo_FIMO_wo_shapes['methyl_minus'].append(np.nan)
			
			AUROC_from_methyl_wm_FIMO['methyl_plus'].append(TF[1][11])
			AUROC_from_methyl_wm_FIMO['non_sensitive'].append(np.nan)
			AUROC_from_methyl_wm_FIMO['methyl_minus'].append(np.nan)
			
			AUROC_from_methyl_wo_FIMO['methyl_plus'].append(TF[1][12])
			AUROC_from_methyl_wo_FIMO['non_sensitive'].append(np.nan)
			AUROC_from_methyl_wo_FIMO['methyl_minus'].append(np.nan)

		
		elif TF[1][0] <= -1.0:
			AUPRC_wm_FIMO_wm_shapes['methyl_minus'].append(TF[1][1])
			AUPRC_wm_FIMO_wm_shapes['non_sensitive'].append(np.nan)
			AUPRC_wm_FIMO_wm_shapes['methyl_plus'].append(np.nan)
			
			AUPRC_wo_FIMO_wm_shapes['methyl_minus'].append(TF[1][2])
			AUPRC_wo_FIMO_wm_shapes['non_sensitive'].append(np.nan)
			AUPRC_wo_FIMO_wm_shapes['methyl_plus'].append(np.nan)			
		
			AUPRC_wm_FIMO_wo_shapes['methyl_minus'].append(TF[1][3])
			AUPRC_wm_FIMO_wo_shapes['non_sensitive'].append(np.nan)
			AUPRC_wm_FIMO_wo_shapes['methyl_plus'].append(np.nan)
			
			AUPRC_wo_FIMO_wo_shapes['methyl_minus'].append(TF[1][4])
			AUPRC_wo_FIMO_wo_shapes['non_sensitive'].append(np.nan)
			AUPRC_wo_FIMO_wo_shapes['methyl_plus'].append(np.nan)
			
			AUPRC_wm_FIMO['methyl_minus'].append(TF[1][5])
			AUPRC_wm_FIMO['non_sensitive'].append(np.nan)
			AUPRC_wm_FIMO['methyl_plus'].append(np.nan)
			
			AUPRC_wo_FIMO['methyl_minus'].append(TF[1][6])
			AUPRC_wo_FIMO['non_sensitive'].append(np.nan)
			AUPRC_wo_FIMO['methyl_plus'].append(np.nan)

			
			AUROC_wm_FIMO_wm_shapes['methyl_minus'].append(TF[1][7])
			AUROC_wm_FIMO_wm_shapes['non_sensitive'].append(np.nan)
			AUROC_wm_FIMO_wm_shapes['methyl_plus'].append(np.nan)
			
			AUROC_wo_FIMO_wm_shapes['methyl_minus'].append(TF[1][8])
			AUROC_wo_FIMO_wm_shapes['non_sensitive'].append(np.nan)
			AUROC_wo_FIMO_wm_shapes['methyl_plus'].append(np.nan)			
		
			AUROC_wm_FIMO_wo_shapes['methyl_minus'].append(TF[1][9])
			AUROC_wm_FIMO_wo_shapes['non_sensitive'].append(np.nan)
			AUROC_wm_FIMO_wo_shapes['methyl_plus'].append(np.nan)
			
			AUROC_wo_FIMO_wo_shapes['methyl_minus'].append(TF[1][10])
			AUROC_wo_FIMO_wo_shapes['non_sensitive'].append(np.nan)
			AUROC_wo_FIMO_wo_shapes['methyl_plus'].append(np.nan)
			
			AUROC_wm_FIMO['methyl_minus'].append(TF[1][11])
			AUROC_wm_FIMO['non_sensitive'].append(np.nan)
			AUROC_wm_FIMO['methyl_plus'].append(np.nan)
			
			AUROC_wo_FIMO['methyl_minus'].append(TF[1][12])
			AUROC_wo_FIMO['non_sensitive'].append(np.nan)
			AUROC_wo_FIMO['methyl_plus'].append(np.nan)
				
			
		else:
			AUPRC_wm_FIMO_wm_shapes['non_sensitive'].append(TF[1][1])
			AUPRC_wm_FIMO_wm_shapes['methyl_plus'].append(np.nan)
			AUPRC_wm_FIMO_wm_shapes['methyl_minus'].append(np.nan)
			
			AUPRC_wo_FIMO_wm_shapes['non_sensitive'].append(TF[1][2])
			AUPRC_wo_FIMO_wm_shapes['methyl_plus'].append(np.nan)
			AUPRC_wo_FIMO_wm_shapes['methyl_minus'].append(np.nan)			
		
			AUPRC_wm_FIMO_wo_shapes['non_sensitive'].append(TF[1][3])
			AUPRC_wm_FIMO_wo_shapes['methyl_plus'].append(np.nan)
			AUPRC_wm_FIMO_wo_shapes['methyl_minus'].append(np.nan)
			
			AUPRC_wo_FIMO_wo_shapes['non_sensitive'].append(TF[1][4])
			AUPRC_wo_FIMO_wo_shapes['methyl_plus'].append(np.nan)
			AUPRC_wo_FIMO_wo_shapes['methyl_minus'].append(np.nan)
			
			AUPRC_wm_FIMO['non_sensitive'].append(TF[1][5])
			AUPRC_wm_FIMO['methyl_plus'].append(np.nan)
			AUPRC_wm_FIMO['methyl_minus'].append(np.nan)
			
			AUPRC_wo_FIMO['non_sensitive'].append(TF[1][6])
			AUPRC_wo_FIMO['methyl_plus'].append(np.nan)
			AUPRC_wo_FIMO['methyl_minus'].append(np.nan)

			
			AUROC_wm_FIMO_wm_shapes['non_sensitive'].append(TF[1][7])
			AUROC_wm_FIMO_wm_shapes['methyl_plus'].append(np.nan)
			AUROC_wm_FIMO_wm_shapes['methyl_minus'].append(np.nan)
			
			AUROC_wo_FIMO_wm_shapes['non_sensitive'].append(TF[1][8])
			AUROC_wo_FIMO_wm_shapes['methyl_plus'].append(np.nan)
			AUROC_wo_FIMO_wm_shapes['methyl_minus'].append(np.nan)			
		
			AUROC_wm_FIMO_wo_shapes['non_sensitive'].append(TF[1][9])
			AUROC_wm_FIMO_wo_shapes['methyl_plus'].append(np.nan)
			AUROC_wm_FIMO_wo_shapes['methyl_minus'].append(np.nan)
			
			AUROC_wo_FIMO_wo_shapes['non_sensitive'].append(TF[1][10])
			AUROC_wo_FIMO_wo_shapes['methyl_plus'].append(np.nan)
			AUROC_wo_FIMO_wo_shapes['methyl_minus'].append(np.nan)
			
			AUROC_wm_FIMO['non_sensitive'].append(TF[1][11])
			AUROC_wm_FIMO['methyl_plus'].append(np.nan)
			AUROC_wm_FIMO['methyl_minus'].append(np.nan)
			
			AUROC_wo_FIMO['non_sensitive'].append(TF[1][12])
			AUROC_wo_FIMO['methyl_plus'].append(np.nan)
			AUROC_wo_FIMO['methyl_minus'].append(np.nan)
	
	
	return AUPRC_wm_FIMO_wm_shapes, AUPRC_wo_FIMO_wm_shapes, AUPRC_wm_FIMO_wo_shapes, AUPRC_wo_FIMO_wo_shapes, AUPRC_wm_FIMO, AUPRC_wo_FIMO, AUROC_wm_FIMO_wm_shapes, AUROC_wo_FIMO_wm_shapes, AUROC_wm_FIMO_wo_shapes, AUROC_wo_FIMO_wo_shapes, AUROC_wm_FIMO, AUROC_wo_FIMO
	

def plot_against(ax, X, Y):
	"""
	'ax' is an axes instance (from matplotlib.pyplot), 'X' and 'Y' are dictionaries containing 3 lists: these contain the scores of the TF binding prediction for a given method (using FIMO with methylation or not, using DNAshapedTFBS or not and including methylation data in DNAshapedTFBS or not). The lists correspond to the methylation sensitivity classification of the TFs. Each list is ordered by the methylation sensitivity of the TFs and contains NA values for every TF that doesn't correspond to that list.
	
	This function plots the scores on the axes instance.
	"""
	
	msize = 4
	lw = 0.5
	
	ax.scatter(X['methyl_minus'], Y['methyl_minus'], s=msize, linewidths=lw, marker='o', edgecolors='limegreen', facecolors='none', alpha=1, label='Methyl -')
	ax.scatter(X['non_sensitive'], Y['non_sensitive'], s=msize, linewidths=lw, marker='o', edgecolors='crimson', facecolors='none', alpha=1, label='Non-sensitive')
	ax.scatter(X['methyl_plus'], Y['methyl_plus'], s=msize, linewidths=lw, marker='o', edgecolors='navy', facecolors='none', alpha=1, label='Methyl +')
	
	ax.tick_params(axis='both', which='both', bottom='on', top='off', left='on', right='off', labelbottom='on', labelleft='on')
#	ax.set_xticks(np.arange(0.0, 1.01, 0.1), minor=False)
#	ax.set_xticks(np.arange(0.0, 1.01, 0.05), minor=True)
#	ax.set_yticks(np.arange(0.0, 1.01, 0.05), minor=True)
	ax.set_xlim(0.0, 1.0)
	ax.set_ylim(0.0, 1.0)
	ax.plot([0.0,1], [0.0,1], color='darkorange', linewidth=lw, linestyle='--')
	ax.legend(loc='lower right', fontsize='xx-small', edgecolor='darkgrey')
	ax.set_aspect('equal', adjustable='box-forced', anchor='C')
	


def plot_graph_combined(ax, wm_FIMO_wm_shapes, wo_FIMO_wm_shapes, wm_FIMO_wo_shapes, wo_FIMO_wo_shapes, wm_FIMO, wo_FIMO):
	"""
	'ax' is an axes instance (from matplotlib.pyplot) and the other inputs are dictionaries containing 3 lists each, differentiating the TFs by methylation sensitivity. Each one of these 3 lists contain the scores of TF binding predictions for the given method and methylation data used.
		
	This function just plots the scores on the axes instance
	"""
	
	msize = 4
	lw = 0.5
	
	xticks = [float(i)/(len(wm_FIMO_wm_shapes['methyl_plus'])+1.0) for i in range(1,len(wm_FIMO_wm_shapes['methyl_plus'])+1)]
	
	ax.scatter(xticks, wm_FIMO_wm_shapes['methyl_plus'], s=msize, linewidths=lw, marker='o', edgecolors='navy', facecolors='none', label='Methyl +, FIMO wm, shapes wm')
	ax.scatter(xticks, wo_FIMO_wm_shapes['methyl_plus'], s=msize, linewidths=lw, marker='^', edgecolors='mediumblue', facecolors='none', label='Methyl +, FIMO wo, shapes wm')
	ax.scatter(xticks, wm_FIMO_wo_shapes['methyl_plus'], s=msize, linewidths=lw, marker='X', edgecolors='blue', facecolors='none', label='Methyl +, FIMO wm, shapes wo')
	ax.scatter(xticks, wo_FIMO_wo_shapes['methyl_plus'], s=msize, linewidths=lw, marker='*', edgecolors='royalblue', facecolors='none', label='Methyl +, FIMO wo, shapes wo')
	ax.scatter(xticks, wm_FIMO['methyl_plus'], s=msize, linewidths=lw, marker='s', edgecolors='deepskyblue', facecolors='none', label='Methyl +, FIMO wm')
	ax.scatter(xticks, wo_FIMO['methyl_plus'], s=msize, linewidths=lw, marker='P', edgecolors='lightskyblue', facecolors='none', label='Methyl +, FIMO wo')
	
	ax.scatter(xticks, wm_FIMO_wm_shapes['non_sensitive'], s=msize, linewidths=lw, marker='o', edgecolors='darkred', facecolors='none', label='Non-sensitive, FIMO wm, shapes wm')
	ax.scatter(xticks, wo_FIMO_wm_shapes['non_sensitive'], s=msize, linewidths=lw, marker='^', edgecolors='red', facecolors='none', label='Non-sensitive, FIMO wo, shapes wm')
	ax.scatter(xticks, wm_FIMO_wo_shapes['non_sensitive'], s=msize, linewidths=lw, marker='X', edgecolors='indianred', facecolors='none', label='Non-sensitive, FIMO wm, shapes wo')
	ax.scatter(xticks, wo_FIMO_wo_shapes['non_sensitive'], s=msize, linewidths=lw, marker='*', edgecolors='tomato', facecolors='none', label='Non-sensitive, FIMO wo, shapes wo')
	ax.scatter(xticks, wm_FIMO['non_sensitive'], s=msize, linewidths=lw, marker='s', edgecolors='darkorange', facecolors='none', label='Non-sensitive, FIMO wm')
	ax.scatter(xticks, wo_FIMO['non_sensitive'], s=msize, linewidths=lw, marker='P', edgecolors='sandybrown', facecolors='none', label='Non-sensitive, FIMO wo')
	
	ax.scatter(xticks, wm_FIMO_wm_shapes['methyl_minus'], s=msize, linewidths=lw, marker='o', edgecolors='darkgreen', facecolors='none', label='Methyl -, FIMO wm, shapes wm')
	ax.scatter(xticks, wo_FIMO_wm_shapes['methyl_minus'], s=msize, linewidths=lw, marker='^', edgecolors='forestgreen', facecolors='none', label='Methyl -, FIMO wo, shapes wm')
	ax.scatter(xticks, wm_FIMO_wo_shapes['methyl_minus'], s=msize, linewidths=lw, marker='X', edgecolors='limegreen', facecolors='none', label='Methyl -, FIMO wm, shapes wo')
	ax.scatter(xticks, wo_FIMO_wo_shapes['methyl_minus'], s=msize, linewidths=lw, marker='*', edgecolors='lime', facecolors='none', label='Methyl -, FIMO wo, shapes wo')
	ax.scatter(xticks, wm_FIMO['methyl_minus'], s=msize, linewidths=lw, marker='s', edgecolors='lightgreen', facecolors='none', label='Methyl -, FIMO wm')
	ax.scatter(xticks, wo_FIMO['methyl_minus'], s=msize, linewidths=lw, marker='P', edgecolors='turquoise', facecolors='none', label='Methyl -, FIMO wo')
		
	ax.set_xticklabels(xticks)
	ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
	ax.set_yticks(np.arange(0.0, 1.0, 0.05), minor=True)
	ax.set_xlim(0.0, 1.0)
	ax.set_ylim(0.0, 1.0)
	ax.legend(bbox_to_anchor=(0.995, 0.75), fontsize='xx-small', edgecolor='darkgrey', bbox_transform=plt.gcf().transFigure)
	ax.set_aspect(2, adjustable='box-forced', anchor='W')
	
	

def plot_violin(ax, X, Y):
	"""
	'ax' is an axes instance (from matplotlib.pyplot), 'X' and 'Y' are dictionaries containing 3 lists: these contain the scores of the TF binding prediction for a given method (using FIMO with methylation or not, using DNAshapedTFBS or not and including methylation data in DNAshapedTFBS or not). The lists correspond to the methylation sensitivity classification of the TFs. Each list is ordered by the methylation sensitivity of the TFs and contains NA values for every TF that doesn't correspond to that list.
	
	This function plots the scores on the axes instance.
	"""
	
	msize = 4
	lw = 0.2
	
	X_arrays = {'methyl_minus': np.asarray([val for val in X['methyl_minus'] if not np.isnan(val)]), 'non_sensitive': np.asarray([val for val in X['non_sensitive'] if not np.isnan(val)]), 'methyl_plus': np.asarray([val for val in X['methyl_plus'] if not np.isnan(val)])}
	Y_arrays = {'methyl_minus': np.asarray([val for val in Y['methyl_minus'] if not np.isnan(val)]), 'non_sensitive': np.asarray([val for val in Y['non_sensitive'] if not np.isnan(val)]), 'methyl_plus': np.asarray([val for val in Y['methyl_plus'] if not np.isnan(val)])}
	
	values = [Y_arrays['methyl_minus']-X_arrays['methyl_minus'], Y_arrays['non_sensitive']-X_arrays['non_sensitive'], Y_arrays['methyl_plus']-X_arrays['methyl_plus']]
	
	ax.violinplot(values, [-1, 0, 1], widths=0.5, points=200)
	ax.plot([-1.5,1.5],[0,0], color='darkorange', linewidth=lw, linestyle='--')
	ax.set_yticks(np.arange(-0.2, 0.2, 0.05), minor=True)
	ax.set_xlim(-1.5, 1.5)
	ax.set_ylim(-0.2, 0.2)
	ax.set_aspect(8, adjustable='box-forced', anchor='C')
	


def generate_graphs(AUPRC_file, AUROC_file, methyl_sensitivity_file):
	
	location = os.path.dirname(AUPRC_file)

	scores = parse_scores(AUPRC_file, AUROC_file, methyl_sensitivity_file)

	AUPRC_wm_FIMO_wm_shapes, AUPRC_wo_FIMO_wm_shapes, AUPRC_wm_FIMO_wo_shapes, AUPRC_wo_FIMO_wo_shapes, AUPRC_wm_FIMO, AUPRC_wo_FIMO, AUROC_wm_FIMO_wm_shapes, AUROC_wo_FIMO_wm_shapes, AUROC_wm_FIMO_wo_shapes, AUROC_wo_FIMO_wo_shapes, AUROC_wm_FIMO, AUROC_wo_FIMO = classify_by_sensitivity(scores)
	
	x_labelpad = 5
	labels_fontsize=10
	title_fontsize=15
	title_pos = 0.95
	shift_left = 0.01
	shift_right = 0.03
	
	
	AUPRC1, prc_ax1 = plt.subplots(nrows=1, ncols=2, num=1)
	AUPRC2, prc_ax2 = plt.subplots(nrows=1, ncols=2, num=2)
	AUPRC3, prc_ax3 = plt.subplots(nrows=1, ncols=2, num=3)
	AUPRC4, prc_ax4 = plt.subplots(nrows=1, ncols=2, num=4)
	AUPRC5, prc_ax5 = plt.subplots(nrows=1, ncols=2, num=5)
	
	AUPRC_violin, prc_ax6 = plt.subplots(nrows=1, ncols=2, num=6)
	
	pos_l = prc_ax1[0].get_position()
	pos_r = prc_ax1[1].get_position()
		
	pos_l.x0 -= shift_left
	pos_r.x0 += shift_right
	pos_l.x1 -= shift_left
	pos_r.x1 += shift_right
	
	prc_ax1[0].set_position(pos_l)
	prc_ax1[1].set_position(pos_r)
	prc_ax2[0].set_position(pos_l)
	prc_ax2[1].set_position(pos_r)
	prc_ax3[0].set_position(pos_l)
	prc_ax3[1].set_position(pos_r)
	prc_ax4[0].set_position(pos_l)
	prc_ax4[1].set_position(pos_r)
	prc_ax5[0].set_position(pos_l)
	prc_ax5[1].set_position(pos_r)
	
	AUROC1, roc_ax1 = plt.subplots(nrows=1, ncols=2, num=7)
	AUROC2, roc_ax2 = plt.subplots(nrows=1, ncols=2, num=8)
	AUROC3, roc_ax3 = plt.subplots(nrows=1, ncols=2, num=9)
	AUROC4, roc_ax4 = plt.subplots(nrows=1, ncols=2, num=10)
	AUROC5, roc_ax5 = plt.subplots(nrows=1, ncols=2, num=11)
	
	AUROC_violin, roc_ax6 = plt.subplots(nrows=1, ncols=2, num=12)
	
	roc_ax1[0].set_position(pos_l)
	roc_ax1[1].set_position(pos_r)
	roc_ax2[0].set_position(pos_l)
	roc_ax2[1].set_position(pos_r)
	roc_ax3[0].set_position(pos_l)
	roc_ax3[1].set_position(pos_r)
	roc_ax4[0].set_position(pos_l)
	roc_ax4[1].set_position(pos_r)
	roc_ax5[0].set_position(pos_l)
	roc_ax5[1].set_position(pos_r)
	
	pos_l.x0 += shift_left
	pos_r.x0 += 1.5*shift_right
	pos_l.x1 += shift_left
	pos_r.x1 += 1.5*shift_right
		
	prc_ax6[0].set_position(pos_l)
	prc_ax6[1].set_position(pos_r)
	
	roc_ax6[0].set_position(pos_l)
	roc_ax6[1].set_position(pos_r)
	
	
#AUPRC
	
	plot_violin(prc_ax6[0], AUPRC_wo_FIMO_wo_shapes, AUPRC_wm_FIMO_wm_shapes)
	prc_ax6[0].set_xlabel("Methylation sensitivity", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax6[0].set_ylabel("FIMO + DNA shapes (w/ - w/o)", fontsize=labels_fontsize)
	
	plot_violin(prc_ax6[1], AUPRC_wo_FIMO, AUPRC_wm_FIMO)
	prc_ax6[1].set_xlabel("Methylation sensitivity", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax6[1].set_ylabel("FIMO scores (w/ - w/o)", fontsize=labels_fontsize)
	
	AUPRC_violin.suptitle("Comparison of predictions AUPRC scores\nwith and without methylation data", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUPRC_violin.savefig("{0}/AUPRC_violin.png".format(location), format='png', dpi=600)

	
	
	plot_against(prc_ax1[0], AUPRC_wo_FIMO_wo_shapes, AUPRC_wm_FIMO_wm_shapes)
	prc_ax1[0].set_xlabel("FIMO w/o Shapes w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax1[0].set_ylabel("FIMO w/ Shapes w/", fontsize=labels_fontsize)
	
	plot_against(prc_ax1[1], AUPRC_wo_FIMO, AUPRC_wm_FIMO)
	prc_ax1[1].set_xlabel("FIMO w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax1[1].set_ylabel("FIMO w/", fontsize=labels_fontsize)
	
	AUPRC1.suptitle("Comparison of predictions AUPRC scores\nwith and without methylation data", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUPRC1.savefig("{0}/AUPRC1.png".format(location), format='png', dpi=600)
	
	plot_against(prc_ax2[0], AUPRC_wo_FIMO_wo_shapes, AUPRC_wo_FIMO_wm_shapes)
	prc_ax2[0].set_xlabel("FIMO w/o Shapes w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax2[0].set_ylabel("FIMO w/o Shapes w/", fontsize=labels_fontsize)
		
	plot_against(prc_ax2[1], AUPRC_wm_FIMO_wo_shapes, AUPRC_wm_FIMO_wm_shapes)
	prc_ax2[1].set_xlabel("FIMO w/ Shapes w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax2[1].set_ylabel("FIMO w/ Shapes w/", fontsize=labels_fontsize)
	
	AUPRC2.suptitle("Comparison of predictions AUPRC scores\nshowing the contribution of\nmethylation data in DNA shapes", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUPRC2.savefig("{0}/AUPRC2.png".format(location), format='png', dpi=600)
	
	
	plot_against(prc_ax3[0], AUPRC_wo_FIMO_wm_shapes, AUPRC_wm_FIMO_wm_shapes)
	prc_ax3[0].set_xlabel("FIMO w/o Shapes w/", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax3[0].set_ylabel("FIMO w/ Shapes w/", fontsize=labels_fontsize)
	
	plot_against(prc_ax3[1], AUPRC_wo_FIMO_wo_shapes, AUPRC_wm_FIMO_wo_shapes)
	prc_ax3[1].set_xlabel("FIMO w/o Shapes w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax3[1].set_ylabel("FIMO w/ Shapes w/o", fontsize=labels_fontsize)
	
	AUPRC3.suptitle("Comparison of predictions AUPRC scores\nshowing the contribution of\nmethylation data in motif and FIMO", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUPRC3.savefig("{0}/AUPRC3.png".format(location), format='png', dpi=600)
	
	
	plot_against(prc_ax4[0], AUPRC_wo_FIMO, AUPRC_wo_FIMO_wo_shapes)
	prc_ax4[0].set_xlabel("FIMO w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax4[0].set_ylabel("FIMO w/o Shapes w/o", fontsize=labels_fontsize)
	
	plot_against(prc_ax4[1], AUPRC_wo_FIMO, AUPRC_wo_FIMO_wm_shapes)
	prc_ax4[1].set_xlabel("FIMO w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax4[1].set_ylabel("FIMO w/o Shapes w/", fontsize=labels_fontsize)
	
	AUPRC4.suptitle("Comparison of predictions AUPRC scores\nto FIMO alone\nwithout methylation data", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUPRC4.savefig("{0}/AUPRC4.png".format(location), format='png', dpi=600)
	
	
	plot_against(prc_ax5[0], AUPRC_wm_FIMO, AUPRC_wm_FIMO_wo_shapes)
	prc_ax5[0].set_xlabel("FIMO w/", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax5[0].set_ylabel("FIMO w/ Shapes w/o", fontsize=labels_fontsize)
	
	plot_against(prc_ax5[1], AUPRC_wm_FIMO, AUPRC_wm_FIMO_wm_shapes)
	prc_ax5[1].set_xlabel("FIMO w/", labelpad=x_labelpad, fontsize=labels_fontsize)
	prc_ax5[1].set_ylabel("FIMO w/ Shapes w/", fontsize=labels_fontsize)
	
	AUPRC5.suptitle("Comparison of predictions AUPRC scores\nto FIMO alone\nwith methylation data", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUPRC5.savefig("{0}/AUPRC5.png".format(location), format='png', dpi=600)
	
	
#AUROC	
	
	plot_violin(roc_ax6[0], AUROC_wo_FIMO_wo_shapes, AUROC_wm_FIMO_wm_shapes)
	roc_ax6[0].set_xlabel("Methylation sensitivity", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax6[0].set_ylabel("FIMO + DNA shapes (w/ - w/o)", fontsize=labels_fontsize)
	
	plot_violin(roc_ax6[1], AUROC_wo_FIMO, AUROC_wm_FIMO)
	roc_ax6[1].set_xlabel("Methylation sensitivity", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax6[1].set_ylabel("FIMO scores (w/ - w/o)", fontsize=labels_fontsize)
	
	AUROC_violin.suptitle("Comparison of predictions AUROC scores\nwith and without methylation data", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUROC_violin.savefig("{0}/AUROC_violin.png".format(location), format='png', dpi=600)
	
		
	plot_against(roc_ax1[0], AUROC_wo_FIMO_wo_shapes, AUROC_wm_FIMO_wm_shapes)
	roc_ax1[0].set_xlabel("FIMO w/o Shapes w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax1[0].set_ylabel("FIMO w/ Shapes w/", fontsize=labels_fontsize)
		
	plot_against(roc_ax1[1], AUROC_wo_FIMO, AUROC_wm_FIMO)
	roc_ax1[1].set_xlabel("FIMO w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax1[1].set_ylabel("FIMO w/", fontsize=labels_fontsize)
	
	AUROC1.suptitle("Comparison of predictions AUROC scores\nwith and without methylation data", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUROC1.savefig("{0}/AUROC1.png".format(location), format='png', dpi=600)
	
	plot_against(roc_ax2[0], AUROC_wo_FIMO_wo_shapes, AUROC_wo_FIMO_wm_shapes)
	roc_ax2[0].set_xlabel("FIMO w/o Shapes w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax2[0].set_ylabel("FIMO w/o Shapes w/", fontsize=labels_fontsize)
	
	plot_against(roc_ax2[1], AUROC_wm_FIMO_wo_shapes, AUROC_wm_FIMO_wm_shapes)
	roc_ax2[1].set_xlabel("FIMO w/ Shapes w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax2[1].set_ylabel("FIMO w/ Shapes w/", fontsize=labels_fontsize)
		
	AUROC2.suptitle("Comparison of predictions AUPRC scores\nshowing the contribution of\nmethylation data in DNA shapes", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUROC2.savefig("{0}/AUROC2.png".format(location), format='png', dpi=600)

	
	plot_against(roc_ax3[0], AUROC_wo_FIMO_wm_shapes, AUROC_wm_FIMO_wm_shapes)
	roc_ax3[0].set_xlabel("FIMO w/o Shapes w/", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax3[0].set_ylabel("FIMO w/ Shapes w/", fontsize=labels_fontsize)
	
	plot_against(roc_ax3[1], AUROC_wo_FIMO_wo_shapes, AUROC_wm_FIMO_wo_shapes)
	roc_ax3[1].set_xlabel("FIMO w/o Shapes w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax3[1].set_ylabel("FIMO w/ Shapes w/o", fontsize=labels_fontsize)
	
	AUROC3.suptitle("Comparison of predictions AUPRC scores\nshowing the contribution of\nmethylation data in motif and FIMO", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUROC3.savefig("{0}/AUROC3.png".format(location), format='png', dpi=600)
	
	
	plot_against(roc_ax4[0], AUROC_wo_FIMO, AUROC_wo_FIMO_wo_shapes)
	roc_ax4[0].set_xlabel("FIMO w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax4[0].set_ylabel("FIMO w/o Shapes w/o", fontsize=labels_fontsize)
		
	plot_against(roc_ax4[1], AUROC_wo_FIMO, AUROC_wo_FIMO_wm_shapes)
	roc_ax4[1].set_xlabel("FIMO w/o", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax4[1].set_ylabel("FIMO w/o Shapes w/", fontsize=labels_fontsize)
	
	AUROC4.suptitle("Comparison of predictions AUROC scores\nto FIMO alone\nwithout methylation data", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUROC4.savefig("{0}/AUROC4.png".format(location), format='png', dpi=600)
		
	
	plot_against(roc_ax5[0], AUROC_wm_FIMO, AUROC_wm_FIMO_wo_shapes)
	roc_ax5[0].set_xlabel("FIMO w/", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax5[0].set_ylabel("FIMO w/ Shapes w/o", fontsize=labels_fontsize)
	
	plot_against(roc_ax5[1], AUROC_wm_FIMO, AUROC_wm_FIMO_wm_shapes)
	roc_ax5[1].set_xlabel("FIMO w/", labelpad=x_labelpad, fontsize=labels_fontsize)
	roc_ax5[1].set_ylabel("FIMO w/ Shapes w/", fontsize=labels_fontsize)
		
	AUROC5.suptitle("Comparison of predictions AUROC scores\nto FIMO alone\nwith methylation data", fontweight='bold', y=title_pos, fontsize=title_fontsize)
	AUROC5.savefig("{0}/AUROC5.png".format(location), format='png', dpi=600)
	
	
	
	
	
	AUPRC_combined, ax1 = plt.subplots(nrows=1, ncols=1, num=13)
	
	plot_graph_combined(ax1, AUPRC_wm_FIMO_wm_shapes, AUPRC_wo_FIMO_wm_shapes, AUPRC_wm_FIMO_wo_shapes, AUPRC_wo_FIMO_wo_shapes, AUPRC_wm_FIMO, AUPRC_wo_FIMO)
	
	ax1.set_ylabel("AUPRC")
	ax1.set_xlabel("TFs ordered by methyl-sensitivity", labelpad=x_labelpad)
	
	AUPRC_combined.suptitle("AUPRC scores calculated with\ndifferent methods and methylation states", fontweight='bold', y=0.98)
	AUPRC_combined.savefig("{0}/AUPRC_combined.png".format(location), format='png', dpi=600)
	
	
	AUROC_combined, ax2 = plt.subplots(nrows=1, ncols=1, num=14)
	
	plot_graph_combined(ax2, AUROC_wm_FIMO_wm_shapes, AUROC_wo_FIMO_wm_shapes, AUROC_wm_FIMO_wo_shapes, AUROC_wo_FIMO_wo_shapes, AUROC_wm_FIMO, AUROC_wo_FIMO)
	
	ax2.set_ylabel("AUROC")
	ax2.set_xlabel("TFs ordered by methyl-sensitivity", labelpad=x_labelpad)
	
	AUROC_combined.suptitle("AUROC scores calculated with\ndifferent methods and methylation states", fontweight='bold', y=0.98)
	AUROC_combined.savefig("{0}/AUROC_combined.png".format(location), format='png', dpi=600)
	
		
if __name__ == "__main__":
	
	AUPRC_scores = sys.argv[1]
	AUROC_scores = sys.argv[2]
	methyl_sensitivity = sys.argv[3]
	
	generate_graphs(AUPRC_scores, AUROC_scores, methyl_sensitivity)
