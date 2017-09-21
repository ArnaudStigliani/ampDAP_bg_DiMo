from calculate_scores_PSSM import calculate_PRC_score
from calculate_scores_PSSM import calculate_ROC_score

def plot_PRC_arg_parsing(subparsers):
	""" Calculate AUPRC score and plot Precision Recall curve. """
	help_str = "Calculate AUPRC score and plot Precision Recall curve."
	parser_prc = subparsers.add_parser('plotPRC', help=help_str)
	parser_prc.add_argument('-f', '--family_directory', required=True, dest='family', 
							action='store', type=str, help='TF-family directory')
	parser_prc.add_argument('-n', '--name_directory', required=True, dest='name', 
							action='store', type=str, help='TF-name directory')
	parser_prc.add_argument('-r', '--results_directory', required=True, dest='results_dir', 
							action='store', type=str, help='results directory')
	parser_prc.add_argument('-t', '--DNAshape_dir', required=True, dest='DNAshape_dir', 
							action='store', type=str, help='directory of DNAshapedTFBS predictions')
	parser_prc.add_argument('-q', '--FIMO_dir', required=True, dest='FIMO_dir', 
							action='store', type=str, help='directory of FIMO predictions')
	parser_prc.add_argument('-P', '--PSSM_predictions', required=False, nargs='?', dest='PSSM_pred', 
							action='store', type=str, choices=['True', 'False'], help='If True, plot the PSSM predictions curve on the same chart. Default is True if the argument is not specified', const='True', default='True')  
	parser_prc.set_defaults(func=calculate_PRC_score)


def plot_ROC_arg_parsing(subparsers):
	""" Calculate AUROC score and plot Receiving Operator curve. """
	help_str = "Calculate AUROC score and plot Receiving Operator curve."
	parser_roc = subparsers.add_parser('plotROC', help=help_str)
	parser_roc.add_argument('-f', '--family_directory', required=True, dest='family', 
							action='store', type=str, help='TF-family directory')
	parser_roc.add_argument('-n', '--name_directory', required=True, dest='name', 
							action='store', type=str, help='TF-name directory')
	parser_roc.add_argument('-r', '--results_directory', required=True, dest='results_dir', 
							action='store', type=str, help='results directory')
	parser_roc.add_argument('-t', '--DNAshape_dir', required=True, dest='DNAshape_dir', 
							action='store', type=str, help='directory of DNAshapedTFBS predictions')
	parser_roc.add_argument('-q', '--FIMO_dir', required=True, dest='FIMO_dir', 
							action='store', type=str, help='directory of FIMO predictions')
	parser_roc.add_argument('-P', '--PSSM_predictions', required=False, nargs='?', dest='PSSM_pred', 
							action='store', type=str, choices=['True', 'False'], help='If True, plot the PSSM predictions curve on the same chart. Default is True if the argument is not specified', const='True', default='True')  
	parser_roc.set_defaults(func=calculate_ROC_score)


def arg_parsing():
	""" Parse the subcommand along with its arguments. """

	descr = '''
	Calculate AUPRC or AUROC score(s) and plot corresponding curve(s).
	'''
	import argparse
	parser = argparse.ArgumentParser(
		description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
	subparsers = parser.add_subparsers(help='Calculate scores and plot corresponding curves',
						title='Subcommands', description='Valid subcommands')
	plot_PRC_arg_parsing(subparsers)
	plot_ROC_arg_parsing(subparsers)
	argu = parser.parse_args()
	return argu
	
