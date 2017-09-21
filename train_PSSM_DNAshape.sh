#!/usr/bin/env bash

#training PSSM + DNA-shape classifier

folder=$1
from_methyl=$2
methyl_shapes=$3
methyl_FIMO=$4

if [ $from_methyl = true ]; then
	subfolder=$folder/from_methyl
	
	if [ $methyl_shapes = true ] && [ $methyl_FIMO = true ]; then
		echo "Training DNAshapedTFBS + PSSM classifier for $subfolder/wm_FIMO_wm_shapes"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		train_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py trainPSSM \
		    -i $folder/foreground/train/fasta/T$indx.fa \
		    -I $folder/foreground/train/bed/T$indx.bed \
		    -b $folder/background/train/fasta/Bt$indx.fa \
		    -B $folder/background/train/bed/Bt$indx.bed \
		    -y -m -t \
		    -g $subfolder/wm_FIMO/predictions/FIMO_predictions_fg_train_$indx.txt \
		    -G $subfolder/wm_FIMO/predictions/FIMO_predictions_bg_train_$indx.txt \
		    -o $subfolder/wm_FIMO_wm_shapes/classifiers/DNAshapedPSSM_classifier_$indx \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Training DNAshapedTFBS + PSSM classifier on $subfolder/wm_FIMO_wm_shapes datasets $indx: $train_time"; 
	'	$folder $subfolder -- {}
	
	elif [ $methyl_shapes = true ] && [ $methyl_FIMO = false ]; then
		echo "Training DNAshapedTFBS + PSSM classifier for $subfolder/wo_FIMO_wm_shapes"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		train_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py trainPSSM \
		    -i $folder/foreground/train/fasta/T$indx.fa \
		    -I $folder/foreground/train/bed/T$indx.bed \
		    -b $folder/background/train/fasta/Bt$indx.fa \
		    -B $folder/background/train/bed/Bt$indx.bed \
		    -y -t \
		    -g $subfolder/wo_FIMO/predictions/FIMO_predictions_fg_train_$indx.txt \
		    -G $subfolder/wo_FIMO/predictions/FIMO_predictions_bg_train_$indx.txt \
		    -o $subfolder/wo_FIMO_wm_shapes/classifiers/DNAshapedPSSM_classifier_$indx \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Training DNAshapedTFBS + PSSM classifier on $subfolder/wo_FIMO_wm_shapes datasets $indx: $train_time"; 
	'	$folder $subfolder -- {}
	
	elif [ $methyl_shapes = false ] && [ $methyl_FIMO = true ]; then
		echo "Training DNAshapedTFBS + PSSM classifier for $subfolder/wm_FIMO_wo_shapes"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		train_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py trainPSSM \
		    -i $folder/foreground/train/fasta/T$indx.fa \
		    -I $folder/foreground/train/bed/T$indx.bed \
		    -b $folder/background/train/fasta/Bt$indx.fa \
		    -B $folder/background/train/bed/Bt$indx.bed \
		    -y -m \
		    -g $subfolder/wm_FIMO/predictions/FIMO_predictions_fg_train_$indx.txt \
		    -G $subfolder/wm_FIMO/predictions/FIMO_predictions_bg_train_$indx.txt \
		    -o $subfolder/wm_FIMO_wo_shapes/classifiers/DNAshapedPSSM_classifier_$indx \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Training DNAshapedTFBS + PSSM classifier on $subfolder/wm_FIMO_wo_shapes datasets $indx: $train_time"; 
	'	$folder $subfolder -- {}
	
	else
		echo "Training DNAshapedTFBS + PSSM classifier for $subfolder/wo_FIMO_wo_shapes"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		train_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py trainPSSM \
		    -i $folder/foreground/train/fasta/T$indx.fa \
		    -I $folder/foreground/train/bed/T$indx.bed \
		    -b $folder/background/train/fasta/Bt$indx.fa \
		    -B $folder/background/train/bed/Bt$indx.bed \
		    -y \
		    -g $subfolder/wo_FIMO/predictions/FIMO_predictions_fg_train_$indx.txt \
		    -G $subfolder/wo_FIMO/predictions/FIMO_predictions_bg_train_$indx.txt \
		    -o $subfolder/wo_FIMO_wo_shapes/classifiers/DNAshapedPSSM_classifier_$indx \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Training DNAshapedTFBS + PSSM classifier on $subfolder/wo_FIMO_wo_shapes datasets $indx: $train_time"; 
	'	$folder $subfolder -- {}
	fi
else
	subfolder=$folder/from_non_methyl

	if [ $methyl_shapes = true ] && [ $methyl_FIMO = true ]; then
		echo "Training DNAshapedTFBS + PSSM classifier for $subfolder/wm_FIMO_wm_shapes"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		train_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py trainPSSM \
		    -i $folder/foreground/train/fasta/T$indx.fa \
		    -I $folder/foreground/train/bed/T$indx.bed \
		    -b $folder/background/train/fasta/Bt$indx.fa \
		    -B $folder/background/train/bed/Bt$indx.bed \
		    -m -t \
		    -g $subfolder/wm_FIMO/predictions/FIMO_predictions_fg_train_$indx.txt \
		    -G $subfolder/wm_FIMO/predictions/FIMO_predictions_bg_train_$indx.txt \
		    -o $subfolder/wm_FIMO_wm_shapes/classifiers/DNAshapedPSSM_classifier_$indx \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Training DNAshapedTFBS + PSSM classifier on $subfolder/wm_FIMO_wm_shapes datasets $indx: $train_time"; 
	'	$folder $subfolder -- {}
	
	elif [ $methyl_shapes = true ] && [ $methyl_FIMO = false ]; then
		echo "Training DNAshapedTFBS + PSSM classifier for $subfolder/wo_FIMO_wm_shapes"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		train_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py trainPSSM \
		    -i $folder/foreground/train/fasta/T$indx.fa \
		    -I $folder/foreground/train/bed/T$indx.bed \
		    -b $folder/background/train/fasta/Bt$indx.fa \
		    -B $folder/background/train/bed/Bt$indx.bed \
		    -t \
		    -g $subfolder/wo_FIMO/predictions/FIMO_predictions_fg_train_$indx.txt \
		    -G $subfolder/wo_FIMO/predictions/FIMO_predictions_bg_train_$indx.txt \
		    -o $subfolder/wo_FIMO_wm_shapes/classifiers/DNAshapedPSSM_classifier_$indx \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Training DNAshapedTFBS + PSSM classifier on $subfolder/wo_FIMO_wm_shapes datasets $indx: $train_time"; 
	'	$folder $subfolder -- {}
	
	elif [ $methyl_shapes = false ] && [ $methyl_FIMO = true ]; then
		echo "Training DNAshapedTFBS + PSSM classifier for $subfolder/wm_FIMO_wo_shapes"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		train_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py trainPSSM \
		    -i $folder/foreground/train/fasta/T$indx.fa \
		    -I $folder/foreground/train/bed/T$indx.bed \
		    -b $folder/background/train/fasta/Bt$indx.fa \
		    -B $folder/background/train/bed/Bt$indx.bed \
		    -m \
		    -g $subfolder/wm_FIMO/predictions/FIMO_predictions_fg_train_$indx.txt \
		    -G $subfolder/wm_FIMO/predictions/FIMO_predictions_bg_train_$indx.txt \
		    -o $subfolder/wm_FIMO_wo_shapes/classifiers/DNAshapedPSSM_classifier_$indx \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Training DNAshapedTFBS + PSSM classifier on $subfolder/wm_FIMO_wo_shapes datasets $indx: $train_time"; 
	'	$folder $subfolder -- {}
	
	else
		echo "Training DNAshapedTFBS + PSSM classifier for $subfolder/wo_FIMO_wo_shapes"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		train_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py trainPSSM \
		    -i $folder/foreground/train/fasta/T$indx.fa \
		    -I $folder/foreground/train/bed/T$indx.bed \
		    -b $folder/background/train/fasta/Bt$indx.fa \
		    -B $folder/background/train/bed/Bt$indx.bed \
		    -g $subfolder/wo_FIMO/predictions/FIMO_predictions_fg_train_$indx.txt \
		    -G $subfolder/wo_FIMO/predictions/FIMO_predictions_bg_train_$indx.txt \
		    -o $subfolder/wo_FIMO_wo_shapes/classifiers/DNAshapedPSSM_classifier_$indx \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Training DNAshapedTFBS + PSSM classifier on $subfolder/wo_FIMO_wo_shapes datasets $indx: $train_time"; 
	'	$folder $subfolder -- {}
	fi
fi
