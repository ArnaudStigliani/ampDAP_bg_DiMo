#!/usr/bin/env bash

#launching the DNA-shaped TFBS tool on DAP-seq data

#Define location of the scripts and results directory

../scripts
../results/DAP_DNAshape
../results

read -r -d '' HELP <<EOF
usage: 
program	[-h] -f TF_family -n TF_name
			-p peaks_fg.narrowPeak | -b peaks_bg.narrowPeak
			-g genome
			[-c methyl_genome] [-a alphabet] [-y] [-q] [-t]
			-l chr_lengths [-r results_dir] 
			[-m meme_motif] [-d]


OPTIONS:
	-f	Transcription factor family
	-n	Transcription factor name
	-p	Original foreground peaks file in narrowPeak format (from DAP-seq library)
	-b	Original background peaks file in narrowPeak format (from ampDAP-seq library)
	-g	File containing the whole genome of the organism you want to work with (fasta format)
	-c	File containing the whole methylated genome of the organism you want to work with (fasta format). 
		If not provided, methylation will not be used at all.
	-a	A customized alphabet (e.g. for methylation data) (optional, used with -c and -q and/or -t options)
	-y	Specifies if you want to perform motif discovery on methylated sequences and derive non-methylated motif 
		from it or the other way around. If -y is specified, the motif is retrieved from methylated sequences first,
		if not specified in the other way.
		(optional, requires -c and -a options)
	-q	Specifies if you want to include methylation data for hits discovery, i.e. for motif and FIMO 
		(optional, requires -c and -a options)
	-t	Specifies if you want to include methylation data for DNA shapes 
		(optional, requires -c and -a options)
	-l	File containing the chromosomes' lengths
	-r	(optional) Results directory you want to create the files in 
		(only in case you have to change your usual results directory)
	-m	(optional) Meme file containing the motif. It goes with -d option
	-d	(optional) If present, it means that datasets, files and other directories created in
		"generate_data_files.sh" already exist and they will not be generated again 
		In that case, be sure that everything is present and provide a meme motif (-m) 
	-h	Display this help

To work properly the program needs at least options -f -n -g -l and option -p (or -b and -m). If you want to include methylation data you also have to use options -c, -a and -q and/or -t and/or -y. 
Exiting...
EOF

from_methyl=false
methylation_FIMO=false
methylation_shapes=false
alphabet=false
methyl_genome=false

OPTIND=1

while getopts f:n:p:b:m:g:c:a:yqtl:r:dh flag; do
	case "${flag}"
	in
		f) TF_family=$OPTARG;;
		n) TF_name=$OPTARG;;
		p) peaks_fg=$OPTARG;;
		b) peaks_bg=$OPTARG;;
		m) meme_motif=$OPTARG;;
		g) genome=$OPTARG;;
		c) methyl_genome=$OPTARG;;
		a) alphabet=$OPTARG;;
		y) from_methyl=true;;
		q) methylation_FIMO=true;;
		t) methylation_shapes=true;;
		l) chr_lengths=$OPTARG;;
		r) results_dir=$OPTARG;;
		d) datasets=true;;
		h) echo "$HELP"; exit 1;;
		\?) echo "Wrong usage of the parameters"; echo "$HELP"; exit 1;;
		*) echo "Wrong usage of the parameters"; echo "$HELP"; exit 1;;
	esac
done

if [ -z "$TF_family" ]; then
	echo "You did not provide any TF family."; echo "$HELP"; exit 1
elif [ -z "$TF_name" ]; then
	echo "You did not provide any TF name."; echo "$HELP"; exit 1
elif [ -z "$peaks_fg" ] || [ -z "$peaks_bg" ]; then
	echo "You did not provide the required peaks files for foreground and background (.narrowPeak)."; echo "$HELP"; exit 1
elif [ -z "$genome" ]; then
	echo "You did not provide any reference genome."; echo "$HELP"; exit 1
elif [ -z "$chr_lengths" ]; then
	echo "You did not provide any file referencing the chromosomes' lengths of the corresponding genome."; echo "$HELP"; exit 1
fi

if [ -n "$results_dir" ]; then
	results_data=$results_dir
fi


if [ $datasets ]; then
	if [ -z "$meme_motif" ]; then
		echo "The -d option (datasets already existing) should go with a meme motif (-m option, .txt file)";
		echo "$HELP"; exit 1
	fi
fi


if ([ $from_methyl = true ] || [ $methylation_FIMO = true ] || [ $methylation_shapes = true ]) && ([ $alphabet = false ] || [ $methyl_genome = false ]); then
	echo "To consider methylation data as you asked (-y, -q or -t option) you should provide the methylated genome (-c) and the customized alphabet for methylation (-a)."; echo "$HELP"; exit 1
elif ([ $from_methyl = false ] && [ $methylation_FIMO = false ] && [ $methylation_shapes = false ]) && ([ $alphabet != false ] || [ $methyl_genome != false ]); then
	echo "You provided a methylated genome (-c) or a customized alphabet (-a) for methylation (or both) but did not ask for methylation to be taken into consideration (-y, -q and or -t)."
	echo "$HELP"
	echo -n "Do you want to continue like this anyway ? [y/n] : "
	read answer
	case $answer in
		[Yy]* ) :;;
		[Nn]* ) echo "Exiting..."; exit;;
		* ) echo "Your answer is not among the expected ones. Exiting..."; exit;;
	esac
fi


if [ ! -e $results_shared/commands.txt ]; then 
	touch $results_shared/commands.txt
fi

read -r -d '' COMMANDS <<EOF

"Program launched on `date` from `pwd`

TF family: $TF_family
TF name: $TF_name
Peaks foreground file: $peaks_fg
Peaks background file: $peaks_bg
Genome file: $genome
Genome with methylation: $methyl_genome
Sequences alphabet for methylation: $alphabet
From methylated motif to non-methylated: $from_methyl
Methylation for hits discovery : $methylation_FIMO
Methylation for DNA shapes : $methylation_shapes
File containing the chromosomes' lengths: $chr_lengths
Results are found in: $results_data
"
EOF

echo "$COMMANDS" >> $results_shared/commands.txt

bash $path_scripts/generate_directories.sh -f $TF_family -n $TF_name -r $results_data

folder=$results_data/$TF_family/$TF_name


### Motif discovery and creation of datasets from peaks files (if not alredy existing)
if [ ! $datasets ]; then
	echo "Creating files"
	bash $path_scripts/generate_data_files.sh -f $TF_family -n $TF_name -p $peaks_fg -b $peaks_bg -g $genome -y $from_methyl -q $methylation_FIMO -c $methyl_genome -a $alphabet -l $chr_lengths -r $results_data

	if ([ $from_methyl = true ] && [ $methylation_FIMO = true ]); then
		meme_motif=$results_data/$TF_family/$TF_name/from_methyl/wm_meme/meme_out/meme_mini.txt
		motif_bg=$results_data/$TF_family/$TF_name/from_methyl/wm_meme/background
		meme_alphabet=$results_data/$TF_family/$TF_name/from_methyl/wm_meme/alphabet
	elif ([ $from_methyl = true ] && [ $methylation_FIMO = false ]); then
		meme_motif=$results_data/$TF_family/$TF_name/from_methyl/wo_meme/meme_out/meme_mini.txt
		motif_bg=$results_data/$TF_family/$TF_name/from_non_methyl/wo_meme/background
		meme_alphabet=false
	elif ([ $from_methyl = false ] && [ $methylation_FIMO = true ]); then
		meme_motif=$results_data/$TF_family/$TF_name/from_non_methyl/wm_meme/meme_out/meme_mini.txt
		motif_bg=$results_data/$TF_family/$TF_name/from_methyl/wm_meme/background
		meme_alphabet=$results_data/$TF_family/$TF_name/from_methyl/wm_meme/alphabet
	else
		meme_motif=$results_data/$TF_family/$TF_name/from_non_methyl/wo_meme/meme_out/meme_mini.txt
		motif_bg=$results_data/$TF_family/$TF_name/from_non_methyl/wo_meme/background
		meme_alphabet=false
	fi
fi

tmp=`grep "letter-probability matrix" $meme_motif | awk -F [=] '{print $4}'`
tmp2=${tmp##" "}
nsites=${tmp2%%" E"}




if ([ $from_methyl = true ] && [ $methylation_FIMO = true ]); then
	dimo_folder=$results_data/$TF_family/$TF_name/from_methyl/wm_DiMO
elif ([ $from_methyl = true ] && [ $methylation_FIMO = false ]); then
	dimo_folder=$results_data/$TF_family/$TF_name/from_methyl/wo_DiMO
elif ([ $from_methyl = false ] && [ $methylation_FIMO = true ]); then
	dimo_folder=$results_data/$TF_family/$TF_name/from_non_methyl/wm_DiMO
else
	dimo_folder=$results_data/$TF_family/$TF_name/from_non_methyl/wo_DiMO
fi


meme_matrix=$dimo_folder/meme_matrix.txt
dimo_matrix_in=$dimo_folder/dimo_matrix_in.txt

if [ -e $dimo_matrix_in ]; then
	echo "The initial matrix to compute DiMO optimization, $dimo_matrix_in, already exists, it will be reused"
else
	bash $path_scripts/meme2dimo.sh $meme_motif $meme_matrix $dimo_matrix_in $methylation_FIMO
fi



### Launching DiMO on train sets and then FIMO on train + test ###############################################
### DiMO motif optimization on training datasets #############################################################
### FIMO to find best hits positions for the motif on each sequence ##########################################
ls $folder/foreground/train/bed/T*.bed | xargs -I{} --max-proc=10 bash -c \
../scripts
	
	folder=$0
	dimo_matrix_in=$1
	from_methyl=$2
	methylation_FIMO=$3
	motif_nsites=$4
	motif_bg=$5
	meme_alphabet=$6
	
	tmp=$(basename {}); tmp2=${tmp##"T"}; indx=${tmp2%%".bed"};
		
	fg_train_bed=$folder/foreground/train/bed/T$indx.bed
	bg_train_bed=$folder/background/train/bed/Bt$indx.bed
	fg_test_bed=$folder/foreground/test/bed/F$indx.bed
	bg_test_bed=$folder/background/test/bed/Bf$indx.bed
	
	
	if [ $methylation_FIMO = true ]; then
	
		fg_train_fasta=$folder/foreground/train/methyl_fasta/methyl_T$indx.fa
		bg_train_fasta=$folder/background/train/methyl_fasta/methyl_Bt$indx.fa
		fg_test_fasta=$folder/foreground/test/methyl_fasta/methyl_F$indx.fa
		bg_test_fasta=$folder/background/test/methyl_fasta/methyl_Bf$indx.fa
		
		if [ $from_methyl = true ]; then
			fimo_dir=$folder/from_methyl/wm_FIMO
			dimo_dir=$folder/from_methyl/wm_DiMO
			meme_dir=$folder/from_methyl/wm_meme
		else
			fimo_dir=$folder/from_non_methyl/wm_FIMO
			dimo_dir=$folder/from_non_methyl/wm_DiMO
			meme_dir=$folder/from_non_methyl/wm_meme
		fi
	else
		fg_train_fasta=$folder/foreground/train/fasta/T$indx.fa
		bg_train_fasta=$folder/background/train/fasta/Bt$indx.fa
		fg_test_fasta=$folder/foreground/test/fasta/F$indx.fa
		bg_test_fasta=$folder/background/test/fasta/Bf$indx.fa
		
		if [ $from_methyl = true ]; then
			fimo_dir=$folder/from_methyl/wo_FIMO
			dimo_dir=$folder/from_methyl/wo_DiMO
			meme_dir=$folder/from_methyl/wo_meme
		else
			fimo_dir=$folder/from_non_methyl/wo_FIMO
			dimo_dir=$folder/from_non_methyl/wo_DiMO
			meme_dir=$folder/from_non_methyl/wo_meme
		fi
	fi
	
	fg_train_FIMO_out=$fimo_dir/FIMO_out/FIMO_fg_train_$indx.txt
	bg_train_FIMO_out=$fimo_dir/FIMO_out/FIMO_bg_train_$indx.txt
	fg_test_FIMO_out=$fimo_dir/FIMO_out/FIMO_fg_test_$indx.txt
	bg_test_FIMO_out=$fimo_dir/FIMO_out/FIMO_bg_test_$indx.txt
	
	fg_train_FIMO_pred=$fimo_dir/predictions/FIMO_predictions_fg_train_$indx.txt
	bg_train_FIMO_pred=$fimo_dir/predictions/FIMO_predictions_bg_train_$indx.txt
	fg_test_FIMO_pred=$fimo_dir/predictions/FIMO_predictions_fg_test_$indx.txt
	bg_test_FIMO_pred=$fimo_dir/predictions/FIMO_predictions_bg_test_$indx.txt
	
	dimo_out_name=$dimo_dir/DiMO_$indx
	dimo_matrix_out="$dimo_out_name"_END.pfm
	matrix_out=$dimo_dir/dimo_matrix_$indx.txt
	new_meme_motif=$meme_dir/optimized_motif_$indx.txt
	new_motif_logo=$meme_dir/optimized_motif_logo_$indx.png
	
	if [ -z "$(ls -A $fimo_dir/predictions/)" ]; then  
		
		Rscript $path_scripts/launch_DiMO.R $methylation_FIMO $fg_train_fasta $bg_train_fasta $dimo_matrix_in $dimo_out_name
		
		sed "1d" $dimo_matrix_out | cut -d" " -f3- > $matrix_out
		
		if [ $meme_alphabet = false ]; then
			cat $matrix_out | matrix2meme -dna -numseqs $motif_nsites -bg $motif_bg > $new_meme_motif
		else
			cat $matrix_out | matrix2meme -alph $meme_alphabet -numseqs $motif_nsites -bg $motif_bg > $new_meme_motif 
		fi
		ceqlogo -i1 $new_meme_motif -o $new_motif_logo -f PNG
	
		
		bash $path_scripts/predict_fimo.sh $new_meme_motif $fg_train_bed $fg_train_fasta $fg_train_FIMO_out
		python $path_scripts/reshape_fimo.py $fg_train_FIMO_out $fg_train_FIMO_pred
		
		bash $path_scripts/predict_fimo.sh $new_meme_motif $bg_train_bed $bg_train_fasta $bg_train_FIMO_out
		python $path_scripts/reshape_fimo.py $bg_train_FIMO_out	$bg_train_FIMO_pred
		
		bash $path_scripts/predict_fimo.sh $new_meme_motif $fg_test_bed $fg_test_fasta $fg_test_FIMO_out
		python $path_scripts/reshape_fimo.py $fg_test_FIMO_out $fg_test_FIMO_pred
		
		bash $path_scripts/predict_fimo.sh $new_meme_motif $bg_test_bed $bg_test_fasta $bg_test_FIMO_out
		python $path_scripts/reshape_fimo.py $bg_test_FIMO_out $bg_test_FIMO_pred
		
	else
		echo "FIMO predictions have already been calculated in $fimo_dir"
	fi
						
' $folder $dimo_matrix_in $from_methyl $methylation_FIMO $nsites $motif_bg $meme_alphabet -- {} 




bash $path_scripts/train_PSSM_DNAshape.sh $folder $from_methyl $methylation_shapes $methylation_FIMO
bash $path_scripts/apply_PSSM_DNAshape_fg.sh $folder $from_methyl $methylation_shapes $methylation_FIMO
bash $path_scripts/apply_PSSM_DNAshape_bg.sh $folder $from_methyl $methylation_shapes $methylation_FIMO

if [ $from_methyl = true ] && [ $methylation_shapes = true ] && [ $methylation_FIMO = true ]; then
	python $path_scripts/calculate_scores_PSSM.py plotPRC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_methyl/wm_FIMO_wm_shapes -q $folder/from_methyl/wm_FIMO --PSSM_predictions='True'
	python $path_scripts/calculate_scores_PSSM.py plotROC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_methyl/wm_FIMO_wm_shapes -q $folder/from_methyl/wm_FIMO --PSSM_predictions='True'

elif [ $from_methyl = true ] && [ $methylation_shapes = true ] && [ $methylation_FIMO = false ]; then
	python $path_scripts/calculate_scores_PSSM.py plotPRC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_methyl/wo_FIMO_wm_shapes -q $folder/from_methyl/wo_FIMO --PSSM_predictions='True'
	python $path_scripts/calculate_scores_PSSM.py plotROC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_methyl/wo_FIMO_wm_shapes -q $folder/from_methyl/wo_FIMO --PSSM_predictions='True'

elif [ $from_methyl = true ] && [ $methylation_shapes = false ] && [ $methylation_FIMO = true ]; then
	python $path_scripts/calculate_scores_PSSM.py plotPRC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_methyl/wm_FIMO_wo_shapes -q $folder/from_methyl/wm_FIMO --PSSM_predictions='True'
	python $path_scripts/calculate_scores_PSSM.py plotROC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_methyl/wm_FIMO_wo_shapes -q $folder/from_methyl/wm_FIMO --PSSM_predictions='True'

elif [ $from_methyl = true ] && [ $methylation_shapes = false ] && [ $methylation_FIMO = false ]; then
	python $path_scripts/calculate_scores_PSSM.py plotPRC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_methyl/wo_FIMO_wo_shapes -q $folder/from_methyl/wo_FIMO --PSSM_predictions='True'
	python $path_scripts/calculate_scores_PSSM.py plotROC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_methyl/wo_FIMO_wo_shapes -q $folder/from_methyl/wo_FIMO --PSSM_predictions='True'

elif [ $from_methyl = false ] && [ $methylation_shapes = true ] && [ $methylation_FIMO = true ]; then
	python $path_scripts/calculate_scores_PSSM.py plotPRC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_non_methyl/wm_FIMO_wm_shapes -q $folder/from_non_methyl/wm_FIMO --PSSM_predictions='True'
	python $path_scripts/calculate_scores_PSSM.py plotROC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_non_methyl/wm_FIMO_wm_shapes -q $folder/from_non_methyl/wm_FIMO --PSSM_predictions='True'

elif [ $from_methyl = false ] && [ $methylation_shapes = true ] && [ $methylation_FIMO = false ]; then
	python $path_scripts/calculate_scores_PSSM.py plotPRC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_non_methyl/wo_FIMO_wm_shapes -q $folder/from_non_methyl/wo_FIMO --PSSM_predictions='True'
	python $path_scripts/calculate_scores_PSSM.py plotROC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_non_methyl/wo_FIMO_wm_shapes -q $folder/from_non_methyl/wo_FIMO --PSSM_predictions='True'

elif [ $from_methyl = false ] && [ $methylation_shapes = false ] && [ $methylation_FIMO = true ]; then
	python $path_scripts/calculate_scores_PSSM.py plotPRC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_non_methyl/wm_FIMO_wo_shapes -q $folder/from_non_methyl/wm_FIMO --PSSM_predictions='True'
	python $path_scripts/calculate_scores_PSSM.py plotROC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_non_methyl/wm_FIMO_wo_shapes -q $folder/from_non_methyl/wm_FIMO --PSSM_predictions='True'

else
	python $path_scripts/calculate_scores_PSSM.py plotPRC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_non_methyl/wm_FIMO_wo_shapes -q $folder/from_non_methyl/wm_FIMO --PSSM_predictions='True'
	python $path_scripts/calculate_scores_PSSM.py plotROC \
		-f $TF_family -n $TF_name -r $results_data -t $folder/from_non_methyl/wm_FIMO_wo_shapes -q $folder/from_non_methyl/wm_FIMO --PSSM_predictions='True'
fi
