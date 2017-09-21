#!/usr/bin/env python2.7

### This script calculates the AUROC and AUPRC scores for a particular TF, depending on the methylation states given. It outputs the corresponding graph in the directory of the main files (predictions of the DNAshapedTFBS algorithm) and it prints out the scores in a general file (appening to it if already existing) in the global results directory.
### See the argparsing file for more specific information on how to use it


import os
from sklearn import metrics as skm
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from argparsing import *


def retrieve_predictions(fg_predictions_files, bg_predictions_files, out_pred_file):
	predictions_list = []

	for i in range(0, len(fg_predictions_files)):
		with open(fg_predictions_files[i]) as fg_pred:
			try:
				first_line_fg = fg_pred.readline().rstrip()
				if first_line_fg[-5:] != "proba" and first_line_fg[-5:] != "score":
					raise Exception
			except Exception as e:
				print "There is an unexpected issue with the file {0}. Please check what could be wrong. The first line tells '{1}'".format(os.path.basename(fg_pred.name), first_line_fg)
			else:
				for line in fg_pred:
					predictions_list.append([float(line.split()[5]), '1', line.split()[4]])
		
		with open(bg_predictions_files[i]) as bg_pred:
			try:
				first_line_bg = bg_pred.readline().rstrip()
				if first_line_bg[-5:] != "proba" and first_line_bg[-5:] != "score":
					raise Exception
			except Exception as e:
				print "There is an unexpected issue with the file {0}. Please check what could be wrong. The first line tells '{1}'".format(os.path.basename(bg_pred.name), first_line_bg)
			else:
				for line in bg_pred:
					predictions_list.append([float(line.split()[5]), '0', line.split()[4]])

#The first line of the files must not be taken into consideration (at least in the normal case) but since we checked it with readline() in the try statement it's already out of the list so we don't have to delete it

	file_name = fg_predictions_files[0]
	folder='{0}/'.format(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(file_name)))))))
	file_name = file_name.replace(folder, '')
	file_name = file_name.replace(file_name[-9:], '')
	print "The files {0} seem to be ok".format(file_name)
	
	if not os.path.isfile(out_pred_file):
		with open(out_pred_file, 'w') as pred_file:
			predictions_list_sorted=sorted(predictions_list, key=lambda x: x[0])
			pred_file.write("Proba\tBg/Fg\tSequence\n")	
			for j in range(0, len(predictions_list_sorted)):
				pred_file.write("{0}\t{1}\t{2}\n".format(predictions_list_sorted[j][0], predictions_list_sorted[j][1], predictions_list_sorted[j][2]))

	for i in range(0, len(predictions_list)):
		del predictions_list[i][2]

#creation of a list containing 2 columns numpy arrays with every line from foreground and background. First column is probability, second is 0 if the line is coming from background predictions and 1 for foreground (true labels)
	
	predictions = np.array(predictions_list).astype(np.float)
	return predictions
	

###
#Precision Recall Curve
###

def calculate_PRC_score(argu):
	
	fg_predictions_DNAshape = []
	bg_predictions_DNAshape = []

	for i in range(0, 10):
		fg_predictions_DNAshape.append("{0}/predictions/DNAshapedPSSM_predictions_fg_{1:d}.txt".format(argu.DNAshape_dir, i))
		bg_predictions_DNAshape.append("{0}/predictions/DNAshapedPSSM_predictions_bg_{1:d}.txt".format(argu.DNAshape_dir, i))
	
	pred_DNAshape=retrieve_predictions(fg_predictions_DNAshape, bg_predictions_DNAshape, "{0}/predictions/predictions_DNAshape_PRC.txt".format(argu.DNAshape_dir))

	AUPRC_DNAshape = skm.average_precision_score(pred_DNAshape[:, 1], pred_DNAshape[:, 0])	
	print "AUPRC from FIMO + DNAshape = {0} ({1}/{2}/{3}/{4})".format(AUPRC_DNAshape, argu.family, argu.name, os.path.basename(os.path.dirname(argu.DNAshape_dir)), os.path.basename(argu.DNAshape_dir))

	precision_DNAshape, recall_DNAshape, thresholds_PRC_DNAshape = skm.precision_recall_curve(pred_DNAshape[:, 1], pred_DNAshape[:, 0], pos_label=1)

	plt.figure()
	lw=2
	plt.plot(recall_DNAshape, precision_DNAshape, color='crimson', lw=lw, label='FIMO + DNAshape: AUPRC = {:.3f}'.format(AUPRC_DNAshape))


	if eval(argu.PSSM_pred):
		fg_predictions_PSSM = []
		bg_predictions_PSSM = []
		
		for i in range(0, 10):
			fg_predictions_PSSM.append("{0}/predictions/FIMO_predictions_fg_test_{1:d}.txt".format(argu.FIMO_dir, i))
			bg_predictions_PSSM.append("{0}/predictions/FIMO_predictions_bg_test_{1:d}.txt".format(argu.FIMO_dir, i))
	
		pred_PSSM=retrieve_predictions(fg_predictions_PSSM, bg_predictions_PSSM, "{0}/predictions/predictions_FIMO_PRC.txt".format(argu.FIMO_dir))

		AUPRC_PSSM = skm.average_precision_score(pred_PSSM[:, 1], pred_PSSM[:, 0])	
		print "AUPRC from FIMO only = {0} ({1}/{2}/{3}/{4})".format(AUPRC_PSSM, argu.family, argu.name, os.path.basename(os.path.dirname(argu.FIMO_dir)), os.path.basename(argu.FIMO_dir))

		precision_PSSM, recall_PSSM, thresholds_PRC_PSSM = skm.precision_recall_curve(pred_PSSM[:, 1], pred_PSSM[:, 0], pos_label=1)
																		
		plt.plot(recall_PSSM, precision_PSSM, color='darkorange', lw=lw, label='FIMO: AUPRC = {:.3f}'.format(AUPRC_PSSM))
	
	
	plt.plot([0, 1], [1, 0], color='navy', lw=lw, linestyle='--')
	plt.xlabel("Recall")
	plt.ylabel("Precision")
	plt.title("Precision-Recall Curve\n{0}/{1}/{2}/{3}".format(argu.family, argu.name, os.path.basename(os.path.dirname(argu.DNAshape_dir)), os.path.basename(argu.DNAshape_dir)))
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.legend(loc="lower left")
	plt.gca().set_aspect('equal', adjustable='box-forced', anchor='C')
	plt.savefig("{0}/predictions/PRC_{1}_{2}.png".format(argu.DNAshape_dir, argu.family, argu.name), format='png', dpi=600)

	AUPRC_scores = "{0}/AUPRC_scores.txt".format(argu.results_dir)
	with open(AUPRC_scores, "a") as scores:
		if (os.path.getsize(AUPRC_scores)==0):
			scores.write("TF_family\tTF_name\tInitial_motif\tMethylation_state_&_Method\tScore\n")
		scores.write("{0}\t{1}\t{2}\t{3}\t{4:f}\n".format(argu.family, argu.name, os.path.basename(os.path.dirname(argu.DNAshape_dir)), os.path.basename(argu.DNAshape_dir), AUPRC_DNAshape))
		if eval(argu.PSSM_pred):
			scores.write("{0}\t{1}\t{2}\t{3}\t{4:f}\n".format(argu.family, argu.name, os.path.basename(os.path.dirname(argu.FIMO_dir)), os.path.basename(argu.FIMO_dir), AUPRC_PSSM))
	
###
#Receiving Operator Curve
###

def calculate_ROC_score(argu):
	
	fg_predictions_DNAshape = []
	bg_predictions_DNAshape = []
	
	for i in range(0, 10):
		fg_predictions_DNAshape.append("{0}/predictions/DNAshapedPSSM_predictions_fg_{1:d}.txt".format(argu.DNAshape_dir, i))
		bg_predictions_DNAshape.append("{0}/predictions/DNAshapedPSSM_predictions_bg_{1:d}.txt".format(argu.DNAshape_dir, i))
	
	pred_DNAshape=retrieve_predictions(fg_predictions_DNAshape, bg_predictions_DNAshape, "{0}/predictions/predictions_DNAshape_ROC.txt".format(argu.DNAshape_dir))
	
	AUROC_DNAshape = skm.roc_auc_score(pred_DNAshape[:, 1], pred_DNAshape[:, 0])
	print "AUROC from FIMO + DNAshape = {0} ({1}/{2}/{3}/{4})".format(AUROC_DNAshape, argu.family, argu.name, os.path.basename(os.path.dirname(argu.DNAshape_dir)), os.path.basename(argu.DNAshape_dir))
	
	fpr_DNAshape, tpr_DNAshape, thresholds_ROC_DNAshape = skm.roc_curve(pred_DNAshape[:, 1], pred_DNAshape[:, 0], pos_label=1)
	
	plt.figure()
	lw=2
	plt.plot(fpr_DNAshape, tpr_DNAshape, color='crimson', lw=lw, label='FIMO + DNAshape: AUROC = {:.3f}'.format(AUROC_DNAshape))
	
	
	if eval(argu.PSSM_pred):
		fg_predictions_PSSM = []
		bg_predictions_PSSM = []
	
		for i in range(0, 10):
			fg_predictions_PSSM.append("{0}/predictions/FIMO_predictions_fg_test_{1:d}.txt".format(argu.FIMO_dir, i))
			bg_predictions_PSSM.append("{0}/predictions/FIMO_predictions_bg_test_{1:d}.txt".format(argu.FIMO_dir, i))
	
		pred_PSSM=retrieve_predictions(fg_predictions_PSSM, bg_predictions_PSSM, "{0}/predictions/predictions_FIMO_ROC.txt".format(argu.FIMO_dir))
			
		AUROC_PSSM = skm.roc_auc_score(pred_PSSM[:, 1], pred_PSSM[:, 0])
		print "AUROC from FIMO only = {0} ({1}/{2}/{3}/{4})".format(AUROC_PSSM, argu.family, argu.name, os.path.basename(os.path.dirname(argu.FIMO_dir)), os.path.basename(argu.FIMO_dir))
	
		fpr_PSSM, tpr_PSSM, thresholds_ROC_PSSM = skm.roc_curve(pred_PSSM[:, 1], pred_PSSM[:, 0], pos_label=1)

		plt.plot(fpr_PSSM, tpr_PSSM, color='darkorange', lw=lw, label='FIMO: AUROC = {:.3f}'.format(AUROC_PSSM))
	
	
	plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
	plt.xlabel("False Positive Rate")
	plt.ylabel("True Positive Rate")
	plt.title("Receiving Operator Curve\n{0}/{1}/{2}/{3}".format(argu.family, argu.name, os.path.basename(os.path.dirname(argu.DNAshape_dir)), os.path.basename(argu.DNAshape_dir)))
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.legend(loc="lower right")
	plt.gca().set_aspect('equal', adjustable='box-forced', anchor='C')
	plt.savefig("{0}/predictions/ROC_{1}_{2}.png".format(argu.DNAshape_dir, argu.family, argu.name), format='png', dpi=600)

	AUROC_scores = "{0}/AUROC_scores.txt".format(argu.results_dir)
	with open(AUROC_scores, "a") as scores:
		if (os.path.getsize(AUROC_scores)==0):
			scores.write("TF_family\tTF_name\tInitial_motif\tMethylation_state_&_Method\tScore\n")
		scores.write("{0}\t{1}\t{2}\t{3}\t{4:f}\n".format(argu.family, argu.name, os.path.basename(os.path.dirname(argu.DNAshape_dir)),os.path.basename(argu.DNAshape_dir), AUROC_DNAshape))
		if eval(argu.PSSM_pred):
			scores.write("{0}\t{1}\t{2}\t{3}\t{4:f}\n".format(argu.family, argu.name, os.path.basename(os.path.dirname(argu.FIMO_dir)), os.path.basename(argu.FIMO_dir), AUROC_PSSM))
	
##############################################################################
#                               MAIN
##############################################################################
if __name__ == "__main__":
    arguments = arg_parsing()
    arguments.func(arguments)
