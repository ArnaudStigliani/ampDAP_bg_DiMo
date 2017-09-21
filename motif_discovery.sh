#!/usr/bin/env bash

path_scripts=/storage/mathelierarea/processed/Victor/scripts
results_data=/storage/mathelierarea/processed/Victor/results/DAP_DNAshape
results_shared=/storage/mathelierarea/processed/Victor/results

folder=$1
TF_family=$2
TF_name=$3
from_methyl=$4
genome=$5
alphabet=$6




if [ $from_methyl = true ]; then
	if [ $methylation_FIMO = true ]; then
		if [ -e $folder/from_methyl/wm_meme/meme_out/meme_mini.txt ]; then	
			echo "The methylated motif already exists in $folder/from_methyl/wm_meme/meme_out/ and will be reused from now on."
		else
			meme_start=`date +%s`
			echo "Retrieving motif from $folder/600peaks_wm_'$TF_family'_'$TF_name'.fa"
			meme-chip -oc $folder/from_methyl/wm_meme -meme-nmotifs 1 -dreme-m 1 -alph $alphabet -noecho -spamo-skip -fimo-skip $folder/600peaks_wm_"$TF_family"_"$TF_name".fa
			meme_end=`date +%s`
			meme_time=$((meme_end-meme_start))
			minutes=$((meme_time / 60))
			seconds=$((meme_time % 60))
			echo "Runtime for methylated meme motif discovery = $meme_time s = $minutes min $seconds s for $folder"
			meme2meme $folder/from_methyl/wm_meme/meme_out/meme.txt > $folder/from_methyl/wm_meme/meme_out/meme_mini.txt
		fi
	else
		if [ -e $folder/from_methyl/wo_meme/meme_out/meme_mini.txt ]; then	
			echo "The non-methylated motif already exists in $folder/from_methyl/wo_meme/meme_out/ and will be reused from now on."
		else
			echo "Creating non_methylated motif from the methylated one"
			if [ ! -e $folder/from_methyl/wo_meme ]; then
				mkdir $folder/from_methyl/wo_meme
			fi
			if [ ! -e $folder/from_methyl/wo_meme/meme_out ]; then
				mkdir $folder/from_methyl/wo_meme/meme_out
			fi
#			bash $path_scripts/convert_methyl_motif.sh $folder/from_methyl/wm_meme/meme_out/meme_mini.txt $folder/from_methyl/wm_meme/background $folder/from_methyl/wo_meme/matrix $folder/from_methyl/wo_meme/background $folder/from_methyl/wo_meme/meme_out/meme_mini.txt $folder/from_methyl/wo_meme/meme_out/logo.eps
		fi
	fi
else
	if [ $methylation_FIMO = true ]; then
		if [ -e $folder/from_non_methyl/wm_meme/meme_out/meme_mini.txt ]; then	
			echo "The methylated motif already exists in $folder/from_non_methyl/wm_meme/meme_out/ and will be reused from now on."
		else
#			echo "Creating methylated motif from the non_methylated one"
#			if [ ! -e $folder/from_non_methyl/wo_meme ]; then
#				mkdir $folder/from_non_methyl/wo_meme
#			fi
#			if [ ! -e $folder/from_non_methyl/wo_meme/meme_out ]; then
#				mkdir $folder/from_non_methyl/wo_meme/meme_out
#			fi
#			bash $path_scripts/convert_methyl_motif.sh $folder/from_non_methyl/wo_meme/meme_out/meme_mini.txt $folder/from_non_methyl/wo_meme/background $folder/from_non_methyl/wm_meme/matrix $folder/from_non_methyl/wm_meme/background $folder/from_non_methyl/wm_meme/meme_out/meme_mini.txt $folder/from_non_methyl/wm_meme/meme_out/logo.eps
		fi
	else
		if [ -e $folder/from_non_methyl/wo_meme/meme_out/meme_mini.txt ]; then	
			echo "The non-methylated motif already exists in $folder/from_non_methyl/wo_meme/meme_out/ and will be reused from now on."
		else
			meme_start=`date +%s`
			echo "Retrieving motif from $folder/600peaks_wo_'$TF_family'_'$TF_name'.fa"
			meme-chip -oc $folder/from_non_methyl/wm_meme -meme-nmotifs 1 -dreme-m 1 -dna -noecho -spamo-skip -fimo-skip $folder/600peaks_wo_"$TF_family"_"$TF_name".fa
			meme_end=`date +%s`
			meme_time=$((meme_end-meme_start))
			minutes=$((meme_time / 60))
			seconds=$((meme_time % 60))
			echo "Runtime for methylated meme motif discovery = $meme_time s = $minutes min $seconds s for $folder"
			meme2meme $folder/from_non_methyl/wo_meme/meme_out/meme.txt > $folder/from_non_methyl/wo_meme/meme_out/meme_mini.txt
		fi
	fi
fi
