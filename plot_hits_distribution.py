#!/usr/bin/env python2.7

### This script plots the distribution of the hits positions (for each TF and depending of whether methylation was used or not)
### It takes as input the results directory for the considered TF and creates a graph with 4 histograms representing the distribution of the hits positions for the foreground and background sequences, both with and without methylation. The position retained is the middle of the hit. It also creates a file with all the histogram values (number of hits at each position for the given TF and file) 


import os
import sys
from sklearn import metrics as skm
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from operator import itemgetter


def parse_hits(hits_bed, side_length=50):
	"""
	hits_bed is a list of bed files in which the function will retrieve the hits positions to return them as a list.
	It returns the position of the middle of the hit (float).
	"""
	### return a normalize position (divided by the length of the sequence) ?
	
	hits_pos_strand = []
	hits_neg_strand = []	
	
	for hits_file in hits_bed:
		with open(hits_file, 'r') as hits:
			for line in hits:
				items = line.rstrip().split()
				
				start = int(items[1])											#0-based (BED file)
				end = int(items[2])
				
				seq_start = int(items[3].split(':')[1].split('-')[0])			#Also 0-based
				
				hit_pos = float(start + end - 1)/2 - seq_start - side_length	#We adapt the end to take the middle of the hit	
				
				if items[5] == '+':
					hits_pos_strand.append(hit_pos)
				else:
					hits_neg_strand.append(hit_pos)
	return hits_pos_strand, hits_neg_strand
		

def parse_files(files_directory):
	"""
	Parses the files in the given directory to retrieve the bed hits files and return the list of files names.
	"""
	import glob
	files_names = '{0}/ext_hits_*.bed'.format(files_directory)
	hits_files = glob.glob(files_names)
	
	return hits_files


def plot_graph(ax, hits_pos_strand, hits_neg_strand, side_length = 50):
	"""
	'ax' is an axes instance (from matplotlib.pyplot).
	'hits_pos_strand' and 'hits_neg_strand' are lists of hits positions relatively to the sequence length, 
	on the positive and negative strand respectively.  
	This function just plots the histogram of hits positions on the axes instance.
	"""
	if hits_pos_strand[0].is_integer():
		hist_bins = np.arange(2*side_length + 2)- side_length - 0.5
	else:
		hist_bins = np.arange(2*side_length + 2)- side_length
	
	n, bins_plot, patches = ax.hist([hits_pos_strand, hits_neg_strand], bins = hist_bins, color=['navy','darkorange'], histtype='barstacked', align='mid')
	
	ax.tick_params(axis='both', which='both', bottom='on', top='off', left='on', right='off', labelbottom='on', labelleft='on')
	ax.set_xticks(np.arange(-side_length, side_length+0.1, side_length), minor=False)
	ax.set_xticks(np.arange(-side_length, side_length+0.1, side_length/5), minor=True)
	ax.set_xlim(-side_length, side_length)
	
#	ax.set_aspect('equal', adjustable='box-forced', anchor='C')
	
	n[1] -= n[0]
	
	return n
	
def generate_graphs(TF_directory):
	
	side_length = 50
	
	non_methyl_fg_from_methyl_test_bed = parse_files('{0}/foreground/test/ext_hits_wo_from_methyl'.format(TF_directory))
	non_methyl_fg_from_methyl_test_pos, non_methyl_fg_from_methyl_test_neg = parse_hits(non_methyl_fg_from_methyl_test_bed, side_length)
	
	non_methyl_bg_from_methyl_test_bed = parse_files('{0}/background/test/ext_hits_wo_from_methyl'.format(TF_directory))
	non_methyl_bg_from_methyl_test_pos, non_methyl_bg_from_methyl_test_neg = parse_hits(non_methyl_bg_from_methyl_test_bed, side_length)
	
	methyl_fg_from_methyl_test_bed = parse_files('{0}/foreground/test/ext_hits_wm_from_methyl'.format(TF_directory))
	methyl_fg_from_methyl_test_pos, methyl_fg_from_methyl_test_neg = parse_hits(methyl_fg_from_methyl_test_bed, side_length)
	
	methyl_bg_from_methyl_test_bed = parse_files('{0}/background/test/ext_hits_wm_from_methyl'.format(TF_directory))
	methyl_bg_from_methyl_test_pos, methyl_bg_from_methyl_test_neg = parse_hits(methyl_bg_from_methyl_test_bed, side_length)
	
	
	non_methyl_fg_from_non_methyl_test_bed = parse_files('{0}/foreground/test/ext_hits_wo_from_non_methyl'.format(TF_directory))
	non_methyl_fg_from_non_methyl_test_pos, non_methyl_fg_from_non_methyl_test_neg = parse_hits(non_methyl_fg_from_non_methyl_test_bed, side_length)
	
	non_methyl_bg_from_non_methyl_test_bed = parse_files('{0}/background/test/ext_hits_wo_from_non_methyl'.format(TF_directory))
	non_methyl_bg_from_non_methyl_test_pos, non_methyl_bg_from_non_methyl_test_neg = parse_hits(non_methyl_bg_from_non_methyl_test_bed, side_length)
	
	methyl_fg_from_non_methyl_test_bed = parse_files('{0}/foreground/test/ext_hits_wm_from_non_methyl'.format(TF_directory))
	methyl_fg_from_non_methyl_test_pos, methyl_fg_from_non_methyl_test_neg = parse_hits(methyl_fg_from_non_methyl_test_bed, side_length)
	
	methyl_bg_from_non_methyl_test_bed = parse_files('{0}/background/test/ext_hits_wm_from_non_methyl'.format(TF_directory))
	methyl_bg_from_non_methyl_test_pos, methyl_bg_from_non_methyl_test_neg = parse_hits(methyl_bg_from_non_methyl_test_bed, side_length)
	
	
	fig, axes = plt.subplots(nrows=4, ncols=2, sharex=True, sharey='row', num=1)
	
	non_methyl_fg_from_methyl_test_hist = plot_graph(axes[0][0], non_methyl_fg_from_methyl_test_pos, non_methyl_fg_from_methyl_test_neg, side_length)
	non_methyl_bg_from_methyl_test_hist = plot_graph(axes[1][0], non_methyl_bg_from_methyl_test_pos, non_methyl_bg_from_methyl_test_neg, side_length)
	methyl_fg_from_methyl_test_hist = plot_graph(axes[0][1], methyl_fg_from_methyl_test_pos, methyl_fg_from_methyl_test_neg, side_length)
	methyl_bg_from_methyl_test_hist = plot_graph(axes[1][1], methyl_bg_from_methyl_test_pos, methyl_bg_from_methyl_test_neg, side_length)
	
	non_methyl_fg_from_non_methyl_test_hist = plot_graph(axes[2][0], non_methyl_fg_from_non_methyl_test_pos, non_methyl_fg_from_non_methyl_test_neg, side_length)
	non_methyl_bg_from_non_methyl_test_hist = plot_graph(axes[3][0], non_methyl_bg_from_non_methyl_test_pos, non_methyl_bg_from_non_methyl_test_neg, side_length)
	methyl_fg_from_non_methyl_test_hist = plot_graph(axes[2][1], methyl_fg_from_non_methyl_test_pos, methyl_fg_from_non_methyl_test_neg, side_length)
	methyl_bg_from_non_methyl_test_hist = plot_graph(axes[3][1], methyl_bg_from_non_methyl_test_pos, methyl_bg_from_non_methyl_test_neg, side_length)
	
	
	axes[3][0].set_xlabel("Position in sequence")
	axes[3][1].set_xlabel("Position in sequence")
	axes[0][0].set_ylabel("Foreground hits")
	axes[1][0].set_ylabel("Background hits")
	axes[2][0].set_ylabel("Foreground hits")
	axes[3][0].set_ylabel("Background hits")


	axes[0][0].set_title("W/o methylation")
	axes[0][1].set_title("W/ methylation")
	
	plus_strand = mlines.Line2D([], [], color='navy', lw = 2, label='Strand +')
	minus_strand = mlines.Line2D([], [], color='darkorange', lw = 2, label='Strand -')
	
	plt.legend(handles=[plus_strand, minus_strand], title="Strand of the hit", bbox_to_anchor=(0.995, 0.995), fontsize='x-small', edgecolor='darkgrey', bbox_transform=plt.gcf().transFigure)
	
	
	ax_left_pos = axes[0][0].get_position()
	ax_right_pos = axes[0][1].get_position()
	
	ax_left_pos.y0 -= 0.05
	ax_left_pos.y1 -= 0.05
	ax_right_pos.y0 -= 0.05
	ax_right_pos.y1 -= 0.05
	
	axes[0][0].set_position(ax_left_pos)
	axes[0][1].set_position(ax_right_pos)
	
	fig.suptitle("Hits relative positions on the 101 bp sequences", fontweight='bold', x=0.4, y=0.95)
	fig.savefig("{0}/hits_positions.png".format(TF_directory), format='png', dpi=600)
	
	
	with open('{0}/histogram_values'.format(TF_directory), 'w') as values_file:
		seq_pos = np.arange(-side_length, len(non_methyl_fg_from_methyl_test_hist[0])-side_length)
		
		current_str = '\t'.join(str(i) for i in seq_pos)
		values_file.write('Histogram values by position\t{0}\n'.format(current_str))
		
		current_str = '\t'.join(str(i) for i in non_methyl_fg_from_methyl_test_hist[0])
		values_file.write('Fg test set w/o methylation from methylated motif, Strand +\t{0}\n'.format(current_str))
		current_str = '\t'.join(str(i) for i in non_methyl_fg_from_methyl_test_hist[1])
		values_file.write('Fg test set w/o methylation from methylated motif, Strand -\t{0}\n'.format(current_str))
		
		current_str = '\t'.join(str(i) for i in non_methyl_bg_from_methyl_test_hist[0])
		values_file.write('Bg test set w/o methylation from methylated motif, Strand +\t{0}\n'.format(current_str))
		current_str = '\t'.join(str(i) for i in non_methyl_bg_from_methyl_test_hist[1])
		values_file.write('Bg test set w/o methylation from methylated motif, Strand -\t{0}\n'.format(current_str))
		
		current_str = '\t'.join(str(i) for i in methyl_fg_from_methyl_test_hist[0])
		values_file.write('Fg test set w/ methylation from methylated motif, Strand +\t{0}\n'.format(current_str))
		current_str = '\t'.join(str(i) for i in methyl_fg_from_methyl_test_hist[1])
		values_file.write('Fg test set w/ methylation from methylated motif, Strand -\t{0}\n'.format(current_str))
		
		current_str = '\t'.join(str(i) for i in methyl_bg_from_methyl_test_hist[0])
		values_file.write('Bg test set w/ methylation from methylated motif, Strand +\t{0}\n'.format(current_str))
		current_str = '\t'.join(str(i) for i in methyl_bg_from_methyl_test_hist[1])
		values_file.write('Bg test set w/ methylation from methylated motif, Strand -\t{0}\n'.format(current_str))
		
		
		current_str = '\t'.join(str(i) for i in non_methyl_fg_from_non_methyl_test_hist[0])
		values_file.write('Fg test set w/o methylation from non methylated motif, Strand +\t{0}\n'.format(current_str))
		current_str = '\t'.join(str(i) for i in non_methyl_fg_from_non_methyl_test_hist[1])
		values_file.write('Fg test set w/o methylation from non methylated motif, Strand -\t{0}\n'.format(current_str))
		
		current_str = '\t'.join(str(i) for i in non_methyl_bg_from_non_methyl_test_hist[0])
		values_file.write('Bg test set w/o methylation from non methylated motif, Strand +\t{0}\n'.format(current_str))
		current_str = '\t'.join(str(i) for i in non_methyl_bg_from_non_methyl_test_hist[1])
		values_file.write('Bg test set w/o methylation from non methylated motif, Strand -\t{0}\n'.format(current_str))
		
		current_str = '\t'.join(str(i) for i in methyl_fg_from_non_methyl_test_hist[0])
		values_file.write('Fg test set w/ methylation from non methylated motif, Strand +\t{0}\n'.format(current_str))
		current_str = '\t'.join(str(i) for i in methyl_fg_from_non_methyl_test_hist[1])
		values_file.write('Fg test set w/ methylation from non methylated motif, Strand -\t{0}\n'.format(current_str))
		
		current_str = '\t'.join(str(i) for i in methyl_bg_from_non_methyl_test_hist[0])
		values_file.write('Bg test set w/ methylation from non methylated motif, Strand +\t{0}\n'.format(current_str))
		current_str = '\t'.join(str(i) for i in methyl_bg_from_non_methyl_test_hist[1])
		values_file.write('Bg test set w/ methylation from non methylated motif, Strand -\t{0}\n'.format(current_str))
	
if __name__ == "__main__":
	
	TF_directory = sys.argv[1]
	generate_graphs(TF_directory)
