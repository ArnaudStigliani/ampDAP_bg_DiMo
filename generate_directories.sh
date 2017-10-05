#!/usr/bin/env bash

results_data=../results/DAP_DNAshape

read -r -d '' HELP <<EOF
usage: 
program [-h] -f TF_family -n TF_name -r results_dir		

This script generates all the necessary directories (if they don't already exist) for the datasets, predictions and score results of a particular TF

OPTIONS:
	-f	Transcription factor family
	-n	Transcription factor name
	-r	Results directory you want to create the files in 
	-h	Display this help

Exiting...
EOF

OPTIND=1

while getopts f:n:r:h flag; do
	case "${flag}"
	in
		f) TF_family=$OPTARG;;
		n) TF_name=$OPTARG;;
		r) results_dir=$OPTARG;;
		h) echo "$HELP"; exit 1;;
		\?) echo "Wrong usage of the parameters"; echo "$HELP"; exit 1;;
		*) echo "Wrong usage of the parameters"; echo "$HELP"; exit 1;;
	esac
done

if [ -n "$results_dir" ]; then
	results_data=$results_dir
else
	echo "You should provide a directory to put the results in."
	echo "$HELP"
	exit 1
fi

#Creation of directories 


if [ ! -e $results_data/$TF_family ]; then
	mkdir $results_data/$TF_family
fi
if [ ! -e $results_data/$TF_family/$TF_name ]; then
	mkdir $results_data/$TF_family/$TF_name
fi

#Foreground directories
if [ ! -e $results_data/$TF_family/$TF_name/foreground ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/test ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/test
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/test/bed ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/test/bed
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/test/fasta ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/test/fasta
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/test/methyl_fasta ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/test/methyl_fasta
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/test/ext_hits_wo_from_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/test/ext_hits_wo_from_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/test/ext_hits_wm_from_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/test/ext_hits_wm_from_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/test/ext_hits_wo_from_non_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/test/ext_hits_wo_from_non_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/test/ext_hits_wm_from_non_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/test/ext_hits_wm_from_non_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/train ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/train
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/train/bed ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/train/bed
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/train/fasta ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/train/fasta
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/train/methyl_fasta ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/train/methyl_fasta
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/train/ext_hits_wo_from_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/train/ext_hits_wo_from_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/train/ext_hits_wm_from_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/train/ext_hits_wm_from_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/train/ext_hits_wo_from_non_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/train/ext_hits_wo_from_non_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/foreground/train/ext_hits_wm_from_non_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/foreground/train/ext_hits_wm_from_non_methyl
fi

#Background directories
if [ ! -e $results_data/$TF_family/$TF_name/background ]; then
	mkdir $results_data/$TF_family/$TF_name/background
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/test ]; then
	mkdir $results_data/$TF_family/$TF_name/background/test
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/test/bed ]; then
	mkdir $results_data/$TF_family/$TF_name/background/test/bed
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/test/fasta ]; then
	mkdir $results_data/$TF_family/$TF_name/background/test/fasta
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/test/methyl_fasta ]; then
	mkdir $results_data/$TF_family/$TF_name/background/test/methyl_fasta
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/test/ext_hits_wo_from_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/background/test/ext_hits_wo_from_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/test/ext_hits_wm_from_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/background/test/ext_hits_wm_from_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/test/ext_hits_wo_from_non_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/background/test/ext_hits_wo_from_non_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/test/ext_hits_wm_from_non_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/background/test/ext_hits_wm_from_non_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/train ]; then
	mkdir $results_data/$TF_family/$TF_name/background/train
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/train/bed ]; then
	mkdir $results_data/$TF_family/$TF_name/background/train/bed
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/train/fasta ]; then
	mkdir $results_data/$TF_family/$TF_name/background/train/fasta
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/train/methyl_fasta ]; then
	mkdir $results_data/$TF_family/$TF_name/background/train/methyl_fasta
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/train/ext_hits_wo_from_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/background/train/ext_hits_wo_from_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/train/ext_hits_wm_from_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/background/train/ext_hits_wm_from_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/train/ext_hits_wo_from_non_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/background/train/ext_hits_wo_from_non_methyl
fi
if [ ! -e $results_data/$TF_family/$TF_name/background/train/ext_hits_wm_from_non_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/background/train/ext_hits_wm_from_non_methyl
fi

#Directories when starting from methylated motif
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl
fi

	#Directory for conversion from methylated to unmethylated
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/to_normal ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/to_normal
fi
	#Directories for DiMO optimized motifs
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wm_DiMO ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wm_DiMO
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wo_DiMO ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wo_DiMO
fi
	#FIMO predictions with methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO/FIMO_out ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO/FIMO_out
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO/predictions
fi
	#FIMO predictions without methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO/FIMO_out ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO/FIMO_out
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO/predictions
fi

	#Final predictions FIMO with methylation / shapes with methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wm_shapes ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wm_shapes
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wm_shapes/classifiers ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wm_shapes/classifiers
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wm_shapes/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wm_shapes/predictions
fi

	#Final predictions FIMO without methylation / shapes with methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wm_shapes ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wm_shapes
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wm_shapes/classifiers ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wm_shapes/classifiers
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wm_shapes/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wm_shapes/predictions
fi

	#Final predictions FIMO with methylation / shapes without methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wo_shapes ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wo_shapes
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wo_shapes/classifiers ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wo_shapes/classifiers
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wo_shapes/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wm_FIMO_wo_shapes/predictions
fi

	#Final predictions FIMO without methylation / shapes without methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wo_shapes ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wo_shapes
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wo_shapes/classifiers ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wo_shapes/classifiers
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wo_shapes/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_methyl/wo_FIMO_wo_shapes/predictions
fi



#Directories when starting from non-methylated motif
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl
fi

	#Directory for conversion from normal to methylated
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/to_methyl ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/to_methyl
fi
	#Directories for DiMO optimized motifs
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wm_DiMO ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wm_DiMO
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wo_DiMO ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wo_DiMO
fi
	#FIMO predictions with methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO/FIMO_out ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO/FIMO_out
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO/predictions
fi
	#FIMO predictions without methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO/FIMO_out ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO/FIMO_out
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO/predictions
fi

	#Final predictions FIMO with methylation / shapes with methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wm_shapes ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wm_shapes
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wm_shapes/classifiers ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wm_shapes/classifiers
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wm_shapes/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wm_shapes/predictions
fi

	#Final predictions FIMO without methylation / shapes with methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wm_shapes ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wm_shapes
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wm_shapes/classifiers ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wm_shapes/classifiers
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wm_shapes/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wm_shapes/predictions
fi

	#Final predictions FIMO with methylation / shapes without methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wo_shapes ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wo_shapes
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wo_shapes/classifiers ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wo_shapes/classifiers
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wo_shapes/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wm_FIMO_wo_shapes/predictions
fi

	#Final predictions FIMO without methylation / shapes without methylation
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wo_shapes ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wo_shapes
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wo_shapes/classifiers ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wo_shapes/classifiers
fi
if [ ! -e $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wo_shapes/predictions ]; then
	mkdir $results_data/$TF_family/$TF_name/from_non_methyl/wo_FIMO_wo_shapes/predictions
fi


