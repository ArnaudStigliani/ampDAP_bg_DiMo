#!/usr/bin/env bash

#Applying PSSM + DNA-shape classifier on foreground 

folder=$1
from_methyl=$2
methyl_shapes=$3
methyl_FIMO=$4

if [ $from_methyl = true ]; then
	subfolder=$folder/from_methyl

	if [ $methyl_shapes = true ] && [ $methyl_FIMO = true ]; then
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wm_FIMO_wm_shapes foreground test files"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		test_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py applyPSSM \
		    -i $folder/foreground/test/fasta/F$indx.fa \
		    -I $folder/foreground/test/bed/F$indx.bed \
		    -y -m -t \
		    -g $subfolder/wm_FIMO/predictions/FIMO_predictions_fg_test_$indx.txt \
		    -c $subfolder/wm_FIMO_wm_shapes/classifiers/DNAshapedPSSM_classifier_$indx.pkl \
		    -o $subfolder/wm_FIMO_wm_shapes/predictions/DNAshapedPSSM_predictions_fg_$indx.txt \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n -v 0.0 \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wm_FIMO_wm_shapes foreground test $indx: $test_time"; 
	'	$folder $subfolder -- {}
	
	elif [ $methyl_shapes = true ] && [ $methyl_FIMO = false ]; then
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wo_FIMO_wm_shapes foreground test files"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		test_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py applyPSSM \
		    -i $folder/foreground/test/fasta/F$indx.fa \
		    -I $folder/foreground/test/bed/F$indx.bed \
		    -y -t \
		    -g $subfolder/wo_FIMO/predictions/FIMO_predictions_fg_test_$indx.txt \
		    -c $subfolder/wo_FIMO_wm_shapes/classifiers/DNAshapedPSSM_classifier_$indx.pkl \
		    -o $subfolder/wo_FIMO_wm_shapes/predictions/DNAshapedPSSM_predictions_fg_$indx.txt \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n -v 0.0 \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wo_FIMO_wm_shapes foreground test $indx: $test_time"; 
	'	$folder $subfolder -- {}
	
	elif [ $methyl_shapes = false ] && [ $methyl_FIMO = true ]; then
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wm_FIMO_wo_shapes foreground test files"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		test_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py applyPSSM \
		    -i $folder/foreground/test/fasta/F$indx.fa \
		    -I $folder/foreground/test/bed/F$indx.bed \
		    -y -m \
		    -g $subfolder/wm_FIMO/predictions/FIMO_predictions_fg_test_$indx.txt \
		    -c $subfolder/wm_FIMO_wo_shapes/classifiers/DNAshapedPSSM_classifier_$indx.pkl \
		    -o $subfolder/wm_FIMO_wo_shapes/predictions/DNAshapedPSSM_predictions_fg_$indx.txt \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n -v 0.0 \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wm_FIMO_wo_shapes foreground test $indx: $test_time"; 
	'	$folder $subfolder -- {}
	
	else
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wo_FIMO_wo_shapes foreground test files"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		test_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py applyPSSM \
		    -i $folder/foreground/test/fasta/F$indx.fa \
		    -I $folder/foreground/test/bed/F$indx.bed \
		    -y \
		    -g $subfolder/wo_FIMO/predictions/FIMO_predictions_fg_test_$indx.txt \
		    -c $subfolder/wo_FIMO_wo_shapes/classifiers/DNAshapedPSSM_classifier_$indx.pkl \
		    -o $subfolder/wo_FIMO_wo_shapes/predictions/DNAshapedPSSM_predictions_fg_$indx.txt \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n -v 0.0 \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wo_FIMO_wo_shapes foreground test $indx: $test_time"; 
	'	$folder $subfolder -- {}
	fi
else
	subfolder=$folder/from_non_methyl
	
	if [ $methyl_shapes = true ] && [ $methyl_FIMO = true ]; then
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wm_FIMO_wm_shapes foreground test files"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		test_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py applyPSSM \
		    -i $folder/foreground/test/fasta/F$indx.fa \
		    -I $folder/foreground/test/bed/F$indx.bed \
		    -m -t \
		    -g $subfolder/wm_FIMO/predictions/FIMO_predictions_fg_test_$indx.txt \
		    -c $subfolder/wm_FIMO_wm_shapes/classifiers/DNAshapedPSSM_classifier_$indx.pkl \
		    -o $subfolder/wm_FIMO_wm_shapes/predictions/DNAshapedPSSM_predictions_fg_$indx.txt \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n -v 0.0 \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wm_FIMO_wm_shapes foreground test $indx: $test_time"; 
	'	$folder $subfolder -- {}
	
	elif [ $methyl_shapes = true ] && [ $methyl_FIMO = false ]; then
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wo_FIMO_wm_shapes foreground test files"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		test_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py applyPSSM \
		    -i $folder/foreground/test/fasta/F$indx.fa \
		    -I $folder/foreground/test/bed/F$indx.bed \
		    -t \
		    -g $subfolder/wo_FIMO/predictions/FIMO_predictions_fg_test_$indx.txt \
		    -c $subfolder/wo_FIMO_wm_shapes/classifiers/DNAshapedPSSM_classifier_$indx.pkl \
		    -o $subfolder/wo_FIMO_wm_shapes/predictions/DNAshapedPSSM_predictions_fg_$indx.txt \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n -v 0.0 \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wo_FIMO_wm_shapes foreground test $indx: $test_time"; 
	'	$folder $subfolder -- {}
	
	elif [ $methyl_shapes = false ] && [ $methyl_FIMO = true ]; then
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wm_FIMO_wo_shapes foreground test files"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		test_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py applyPSSM \
		    -i $folder/foreground/test/fasta/F$indx.fa \
		    -I $folder/foreground/test/bed/F$indx.bed \
		    -m \
		    -g $subfolder/wm_FIMO/predictions/FIMO_predictions_fg_test_$indx.txt \
		    -c $subfolder/wm_FIMO_wo_shapes/classifiers/DNAshapedPSSM_classifier_$indx.pkl \
		    -o $subfolder/wm_FIMO_wo_shapes/predictions/DNAshapedPSSM_predictions_fg_$indx.txt \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n -v 0.0 \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wm_FIMO_wo_shapes foreground test $indx: $test_time"; 
	'	$folder $subfolder -- {}
	
	else
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wo_FIMO_wo_shapes foreground test files"
		ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
	'	path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified
	
		folder=$0
		subfolder=$1
		
		tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
		test_time=$((TIMEFORMAT=%3lR; time python2.7 $path_DNAshapedTFBS/DNAshapedTFBS.py applyPSSM \
		    -i $folder/foreground/test/fasta/F$indx.fa \
		    -I $folder/foreground/test/bed/F$indx.bed \
		    -g $subfolder/wo_FIMO/predictions/FIMO_predictions_fg_test_$indx.txt \
		    -c $subfolder/wo_FIMO_wo_shapes/classifiers/DNAshapedPSSM_classifier_$indx.pkl \
		    -o $subfolder/wo_FIMO_wo_shapes/predictions/DNAshapedPSSM_predictions_fg_$indx.txt \
		    -s "HelT" "ProT" "MGW" "Roll" \
		    -n -v 0.0 \
		    1>>"$folder/stdout.txt" 2>>"$folder/stderr.txt";) 2>&1;)
		    
		echo "Applying DNAshapedTFBS + PSSM on $subfolder/wo_FIMO_wo_shapes foreground test $indx: $test_time"; 
	'	$folder $subfolder -- {}
	fi
fi
