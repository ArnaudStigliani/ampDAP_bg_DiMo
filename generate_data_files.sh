#!/usr/bin/env bash

path_scripts=../scripts
path_DNAshapedTFBS=/storage/scratch/victor/DNAshape/src/DNAshapedTFBS_modified	#local
path_BiasAway=/lsc/BiasAway
results_data=../results/DAP_DNAshape
results_shared=../results
data=../data

read -r -d '' HELP <<EOF
usage: 
program [-h] -f TF_family -n TF_name 
		-p peaks_fg.narrowPeak | -b peaks_bg.narrowPeak 
		-g genome 
		[-y bool] [-q methylation_bool]
		[-c methyl_genome] [-a alphabet]
		-l chr_lengths [-r results_dir]

This script generates all the necessary datasets files for foreground and background.

OPTIONS:
	-f	Transcription factor family
	-n	Transcription factor name
	-p	Original peaks foreground file in narrowPeak format (from DAP-seq library)
	-b	Original peaks background file in narrowPeak format (from ampDAP-seq library)
	-g	File containing the whole genome of the organism you want to work with (fasta format)
	-y	Specifies if you want to perform motif discovery on methylated sequences and derive non-methylated motif 
		from it or the other way around. If -y is specified the motif is retrieved from methylated sequences first,
		if not specified in the other way.
		(optional, requires -c and -a options)
	-q	Boolean: specifies whether you want to include methylation data. (if true, requires -c and -a options)
		(optional, default is false)
	-c	File containing the whole methylated genome of the organism you want to work with (fasta format)
		You probably want to use -q and/or -t and -a options with this one
	-a	A customized alphabet (e.g. for methylation data) (optional, used with -c and -q and/or -t options)
	-l	File containing the chromosomes' lengths
	-r	Results directory you want to create the files in 
		(optional, only in case you have to change your usual results directory)
	-h	Display this help

To work properly the program needs options -f -n -g -l, either option -p or option -b (if the bed file already exists)
To include methylation data you need to use option -c, -a and -q
EOF

methylation_FIMO=false
side_length=50

OPTIND=1

while getopts f:n:p:b:g:c:y:q:t:a:l:r:h flag; do
	case "${flag}"
	in
		f) TF_family=$OPTARG;;
		n) TF_name=$OPTARG;;
		p) ori_peaks_fg=$OPTARG;;
		b) peaks_bg=$OPTARG;;
		g) genome=$OPTARG;;
		y) from_methyl=$OPTARG;;
		q) methylation_FIMO=$OPTARG;;
		c) methyl_genome=$OPTARG;;
		a) alphabet=$OPTARG;;
		l) chr_lengths=$OPTARG;;
		r) results_dir=$OPTARG;;
		h) echo "$HELP"; exit 1;;
		\?) echo "Wrong usage of the parameters"; echo "$HELP"; exit 1;;
		*) echo "Wrong usage of the parameters"; echo "$HELP"; exit 1;;
	esac
done

if [ -n "$results_dir" ]; then
	results_data=$results_dir
fi

folder=$results_data/$TF_family/$TF_name
foreground=$folder/foreground
background=$folder/background

#Variable names for files
bed600=$folder/600peaks_"$TF_family"_"$TF_name".bed
fasta600=$folder/600peaks_wo_"$TF_family"_"$TF_name".fa
methyl_fasta600=$folder/600peaks_wm_"$TF_family"_"$TF_name".fa

ori_peaks_fg_amb_bed=$folder/ori_peaks_fg_amb_"$TF_family"_"$TF_name".bed
ori_peaks_fg_seq=$folder/ori_peaks_fg_seq_"$TF_family"_"$TF_name".fa
ori_peaks_fg_bed=$folder/ori_peaks_fg_"$TF_family"_"$TF_name".bed

#Turning narrowPeak file into bed file and sorting it randomly
if [ -n "$ori_peaks_fg" ]; then
	if [ -e $ori_peaks_fg_amb_bed ]; then
		echo "File $(basename $ori_peaks_fg_amb_bed) already exists and will be reused from now on."
	else
		awk ' { $(($3=$2+$NF+'$side_length'+1)); $(($2+=$NF-'$side_length')); print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$5"\t"$6 } ' $ori_peaks_fg | sort -R > $ori_peaks_fg_amb_bed
	fi
fi

#Retrieving normal fasta sequences from bed
if [ -e $ori_peaks_fg_seq ]; then 
	echo "The normal fasta sequences $(basename $ori_peaks_fg_seq) already exist and will not be overwritten."
else
	bedtools getfasta -fi $genome -bed $ori_peaks_fg_amb_bed -fo $folder/ori_peaks_fg_seq_amb_"$TF_family"_"$TF_name".fa

#Removing all sequences that contain any non-ACGT base from foreground (FIMO doesn't handle ambiguous bases)
	python $path_scripts/remove_ambiguous.py $folder/ori_peaks_fg_seq_amb_"$TF_family"_"$TF_name".fa $ori_peaks_fg_seq
fi

#Writing bed file with names corresponding to the fasta file just created and sequences containing ambiguous bases removed
if [ -e $ori_peaks_fg_bed ]; then
	echo "File $(basename $ori_peaks_fg_bed) already exists and will be reused from now on."
else
	grep -e ">" $ori_peaks_fg_seq | cut -c2- | awk -F [:-] -v OFS='\t' '{print $1,$2,$3,$0}' > $ori_peaks_fg_bed
fi







### 600 peaks file creation for motif discovery ####################################################################

if [ -e $bed600 ]; then
	echo "File $(basename $bed600) already exists and will be reused from now on."
else
	head -n 600 $ori_peaks_fg_bed > $bed600
fi
if [ -e $methyl_fasta600 ]; then
	echo "File $(basename $methyl_fasta600) already exists and will be reused from now on."
else
	bedtools getfasta -fi $methyl_genome -bed $bed600 -fo $methyl_fasta600
fi
if [ -e $fasta600 ]; then
	echo "File $(basename $fasta600) already exists and will be reused from now on."
else
	bedtools getfasta -fi $genome -bed $bed600 -fo $fasta600
fi





### MOTIF DISCOVERY ###################################################################################################

if [ $from_methyl = true ]; then
### Starting from methylated motif ####################################################################################
	from_methyl_motif=$folder/from_methyl/wm_meme/meme_out/meme_mini.txt
	methyl_full_motif=$folder/from_methyl/wm_meme/meme_out/meme.txt
	normal_meme_bg=$folder/from_non_methyl/wo_meme/background

	if [ -e $from_methyl_motif ]; then	
		echo "The methylated motif already exists in $folder/from_methyl/wm_meme/meme_out/ and will be reused from now on."
	else
		meme_start=`date +%s`
		echo "Retrieving motif from $methyl_fasta600"
		meme-chip -oc $folder/from_methyl/wm_meme -meme-nmotifs 1 -dreme-m 0 -alph $alphabet -noecho -spamo-skip -fimo-skip $methyl_fasta600
		meme_end=`date +%s`
		meme_time=$((meme_end-meme_start))
		minutes=$((meme_time / 60))
		seconds=$((meme_time % 60))
		echo "Runtime for methylated meme motif discovery = $meme_time s = $minutes min $seconds s for $folder"
		meme2meme $methyl_full_motif > $from_methyl_motif
	fi
		
	if [ $methylation_FIMO = true ]; then
### No need for non-methylated motif if we don't use it with FIMO
		:
	else
### Retrieving unmethylated motif from methylated
		unmethylated_motif=$folder/from_methyl/wo_meme/meme_out/meme_mini.txt
		unmethylated_motif_logo=$folder/from_methyl/wo_meme/meme_out/logo.png
		if [ -e $unmethylated_motif ]; then	
			echo "The non-methylated motif already exists in $folder/from_methyl/wo_meme/meme_out/ and will be reused from now on."
		else
			echo "Creating non_methylated motif from the methylated one"
			if [ ! -e $folder/from_methyl/wo_meme ]; then
				mkdir $folder/from_methyl/wo_meme
			fi
			if [ ! -e $folder/from_methyl/wo_meme/meme_out ]; then
				mkdir $folder/from_methyl/wo_meme/meme_out
			fi

			#Check if background from non-methylated motif exists, else we retrieve non-methylated motif to get its background 
			if [ -e $normal_meme_bg ]; then
				echo "The background from non-methylated motif already exists and will be used"
			else
				from_normal_motif=$folder/from_non_methyl/wo_meme/meme_out/meme_mini.txt
				meme_start=`date +%s`
				echo "Retrieving motif from $fasta600"
				meme-chip -oc $folder/from_non_methyl/wo_meme -meme-nmotifs 1 -dreme-m 0 -dna -noecho -spamo-skip -fimo-skip $fasta600
				meme_end=`date +%s`
				meme_time=$((meme_end-meme_start))
				minutes=$((meme_time / 60))
				seconds=$((meme_time % 60))
				echo "Runtime for non-methylated meme motif discovery = $meme_time s = $minutes min $seconds s for $folder"
				meme2meme $folder/from_non_methyl/wo_meme/meme_out/meme.txt > $from_normal_motif
			fi
		
			#Creation of unmethylated motif from the methylated one
			
			MEME_hits=$folder/from_methyl/to_normal/MEME_hits.txt
			bed_hits=$folder/from_methyl/to_normal/MEME_hits.bed
			fasta_hits=$folder/from_methyl/to_normal/unmethylated_MEME_hits.fa
			sites_file=$folder/from_methyl/to_normal/wm2wo.sites
			
			#Find MEME hits (methylated sequences) to then retrieve non-methylated motif 
			start_line=`grep -n "sites sorted by position p-value" $methyl_full_motif | awk -F [:] '{print $1}'`
			sites_start_line=$((start_line + 4))
			
			tmp=`grep "letter-probability matrix" $methyl_full_motif | awk -F [=] '{print $4}'`
			tmp2=${tmp##" "}
			nsites=${tmp2%%" E"}
			
			sites_end_line=$((sites_start_line + nsites - 1))
			
			head -n $sites_end_line $methyl_full_motif | tail -n $nsites > $MEME_hits
			
			#Turn MEME hits into BED
			awk -v OFS='\t' '{split($1, a, "[:-]"); $((start=a[2]+$3-1)); $((end=a[2]+$3+length($6)-1)); print a[1],start,end,$1,$4,$2}' $MEME_hits > $bed_hits
			
			#Retrieve normal fasta sequences from this BED file to get the non-methylated hits
			bedtools getfasta -s -name -fi $genome -bed $bed_hits -fo $fasta_hits
			
			#Create sites file from fasta (removing headers) and create unmethylated motif from these sites
			sed '/>/d' $fasta_hits > $sites_file
			sites2meme -ext .sites -bg $normal_meme_bg $folder/from_methyl/to_normal > $folder/from_methyl/wo_meme/meme_out/meme.txt
			meme2meme $folder/from_methyl/wo_meme/meme_out/meme.txt > $unmethylated_motif
			ceqlogo -i1 $unmethylated_motif -o $unmethylated_motif_logo -f PNG
		fi
	fi
			
else
### Starting from non-methylated motif ####################################################################################
	from_normal_motif=$folder/from_non_methyl/wo_meme/meme_out/meme_mini.txt
	normal_full_motif=$folder/from_non_methyl/wo_meme/meme_out/meme.txt
	methyl_meme_bg=$folder/from_methyl/wm_meme/background

	if [ -e $from_normal_motif ]; then	
		echo "The non-methylated motif already exists in $folder/from_non_methyl/wo_meme/meme_out/ and will be reused from now on."
	else
		meme_start=`date +%s`
		echo "Retrieving motif from $fasta600"
		meme-chip -oc $folder/from_non_methyl/wo_meme -meme-nmotifs 1 -dreme-m 0 -dna -noecho -spamo-skip -fimo-skip $fasta600
		meme_end=`date +%s`
		meme_time=$((meme_end-meme_start))
		minutes=$((meme_time / 60))
		seconds=$((meme_time % 60))
		echo "Runtime for non-methylated meme motif discovery = $meme_time s = $minutes min $seconds s for $folder"
		meme2meme $normal_full_motif > $from_normal_motif
	fi
		
	if [ $methylation_FIMO = true ]; then
### Retrieving methylated motif from non-methylated
		methylated_motif=$folder/from_non_methyl/wm_meme/meme_out/meme_mini.txt
		methylated_motif_logo=$folder/from_non_methyl/wm_meme/meme_out/logo.png
		if [ -e $methylated_motif ]; then	
			echo "The methylated motif already exists in $folder/from_non_methyl/wm_meme/meme_out/ and will be reused from now on."
		else
			echo "Creating methylated motif from the non-methylated one"
			if [ ! -e $folder/from_non_methyl/wm_meme ]; then
				mkdir $folder/from_non_methyl/wm_meme
			fi
			if [ ! -e $folder/from_non_methyl/wm_meme/meme_out ]; then
				mkdir $folder/from_non_methyl/wm_meme/meme_out
			fi

			#Check if background from original methylated motif exists, else we retrieve methylated motif to get its background 
			if [ -e $methyl_meme_bg ]; then
				echo "The background from non-methylated motif already exists and will be used"
			else
				from_methyl_motif=$folder/from_methyl/wm_meme/meme_out/meme_mini.txt
				meme_start=`date +%s`
				echo "Retrieving motif from $methyl_fasta600"
				meme-chip -oc $folder/from_methyl/wm_meme -meme-nmotifs 1 -dreme-m 0 -alph $alphabet -noecho -spamo-skip -fimo-skip $methyl_fasta600
				meme_end=`date +%s`
				meme_time=$((meme_end-meme_start))
				minutes=$((meme_time / 60))
				seconds=$((meme_time % 60))
				echo "Runtime for methylated meme motif discovery = $meme_time s = $minutes min $seconds s for $folder"
				meme2meme $folder/from_methyl/wm_meme/meme_out/meme.txt > $from_methyl_motif
			fi
		
			#Creation of methylated motif from the non-methylated one
				
			MEME_hits=$folder/from_non_methyl/to_methyl/MEME_hits.txt
			bed_hits=$folder/from_non_methyl/to_methyl/MEME_hits.bed
			k_methyl_fasta_hits=$folder/from_non_methyl/to_methyl/k_methylated_MEME_hits.fa
			methyl_fasta_hits=$folder/from_non_methyl/to_methyl/methylated_MEME_hits.fa
			sites_file=$folder/from_non_methyl/to_methyl/wo2wm.sites
			
			#Find MEME hits (normal sequences) to then retrieve methylated motif 
			start_line=`grep -n "sites sorted by position p-value" $normal_full_motif | awk -F [:] '{print $1}'`
			sites_start_line=$((start_line + 4))
			
			tmp=`grep "letter-probability matrix" $normal_full_motif | awk -F [=] '{print $4}'`
			tmp2=${tmp##" "}
			nsites=${tmp2%%" E"}
			
			sites_end_line=$((sites_start_line + nsites - 1))
			
			head -n $sites_end_line $normal_full_motif | tail -n $nsites > $MEME_hits
			
			#Turn MEME hits into BED
			awk -v OFS='\t' '{split($1, a, "[:-]"); $((start=a[2]+$3-1)); $((end=a[2]+$3+length($6)-1)); print a[1],start,end,$1,$4,$2}' $MEME_hits > $bed_hits
			
			#Retrieve methylated fasta sequences from this BED file to get the methylated hits
			bedtools getfasta -s -name -fi $methyl_genome -bed $bed_hits -fo $k_methyl_fasta_hits
			
			#Adjust methyl_fasta600_hits (problem with complementary bases on minus strand with bedtools)
			python $path_scripts/adjust_comp_methyl.py $k_methyl_fasta_hits $bed_hits $methyl_fasta_hits
			
			#Create sites file from methyl-fasta (removing headers) and create methylated motif from these sites
			sed '/>/d' $methyl_fasta_hits > $sites_file
			sites2meme -ext .sites -alph $alphabet -bg $methyl_meme_bg $folder/from_non_methyl/to_methyl > $folder/from_non_methyl/wm_meme/meme_out/meme.txt
			meme2meme $folder/from_non_methyl/wm_meme/meme_out/meme.txt > $methylated_motif
			ceqlogo -i1 $methylated_motif -o $methylated_motif_logo -f PNG
		fi
	else
### No need for methylated motif if we don't use it with FIMO
		:
	fi
fi



### FOREGROUND from DAP-seq library #################################################################################
### Sequences containing ambiguous bases have already been removed in the previous steps ############################
peaks_fg_bed=$folder/peaks_fg_"$TF_family"_"$TF_name".bed
peaks_fg_seq=$folder/peaks_fg_seq_"$TF_family"_"$TF_name".fa

if [ -e $peaks_fg_bed ]; then
	echo "File $(basename $peaks_fg_bed) already exists and will be reused from now on."
else
	tail -n +601 $ori_peaks_fg_bed > $peaks_fg_bed
fi

#Retrieving normal fasta sequences from bed
if [ -e $peaks_fg_seq ]; then 
	echo "The normal fasta sequences $(basename $peaks_fg_seq) already exist and will not be overwritten."
else
	bedtools getfasta -fi $genome -bed $peaks_fg_bed -fo $peaks_fg_seq
fi


### BACKGROUND from ampDAP-seq library #################################################################################
peaks_bg_all_bed=$folder/peaks_bg_all_"$TF_family"_"$TF_name".bed
peaks_bg_amb_bed=$folder/peaks_bg_amb_"$TF_family"_"$TF_name".bed
peaks_bg_seq=$folder/peaks_bg_seq_"$TF_family"_"$TF_name".fa
peaks_bg_bed=$folder/peaks_bg_"$TF_family"_"$TF_name".bed

#Turning narrowPeak file into bed file and sorting it randomly
if [ -n "$peaks_bg" ]; then
	if [ -e $peaks_bg_all_bed ]; then
		echo "File $(basename $peaks_bg_all_bed) already exists and will be reused from now on."
	else
		awk ' { $(($3=$2+$NF+'$side_length'+1)); $(($2+=$NF-'$side_length')); print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$5"\t"$6 } ' $peaks_bg | sort -R > $peaks_bg_all_bed
	fi
fi

#Retrieving sequences from ampDAP-seq library non-overlapping with DAP-seq to create background
bedtools intersect -v -wa -a $peaks_bg_all_bed -b $ori_peaks_fg_amb_bed > $peaks_bg_amb_bed

#Retrieving normal fasta sequences from bed
if [ -e $peaks_bg_seq ]; then 
	echo "The normal fasta sequences $(basename $peaks_bg_seq) already exist and will not be overwritten."
else
	bedtools getfasta -fi $genome -bed $peaks_bg_amb_bed -fo $folder/peaks_bg_seq_amb_"$TF_family"_"$TF_name".fa

#Removing all sequences that contain any non-ACGT base from background (FIMO doesn't handle ambiguous bases)
	python $path_scripts/remove_ambiguous.py $folder/peaks_bg_seq_amb_"$TF_family"_"$TF_name".fa $peaks_bg_seq
fi

#Writing bed file with names corresponding to the fasta file just created and sequences containing ambiguous bases removed
if [ -e $peaks_bg_bed ]; then
	echo "File $(basename $peaks_bg_bed) already exists and will be reused from now on."
else
	grep -e ">" $peaks_bg_seq | cut -c2- | awk -F [:-] -v OFS='\t' '{print $1,$2,$3,$0}' > $peaks_bg_bed
fi




#Retrieving methyl fasta sequences from bed if necessary
if [ $methylation_FIMO = true ]; then
	peaks_fg_methyl_seq=$folder/peaks_fg_methyl_seq_"$TF_family"_"$TF_name".fa
	peaks_bg_methyl_seq=$folder/peaks_bg_methyl_seq_"$TF_family"_"$TF_name".fa
	
	if [ -e $peaks_fg_methyl_seq ]; then
		echo "The methylated foreground sequences $(basename $peaks_fg_methyl_seq) already exist and will not be overwritten."
	else
		bedtools getfasta -fi $methyl_genome -bed $peaks_fg_bed -fo $peaks_fg_methyl_seq
	fi
	
	if [ -e $peaks_bg_methyl_seq ]; then
		echo "The methylated background sequences $(basename $peaks_bg_methyl_seq) already exist and will not be overwritten."
	else
		bedtools getfasta -fi $methyl_genome -bed $peaks_bg_bed -fo $peaks_bg_methyl_seq
	fi
fi


###############################################################################################################
#Function to create and fill the datasets files (divide the in_file into 10 files with prefix out_file)

chopFile()
{
	in_file=$1
	lpseq=$2			#number of lines per sequence in the file
	out_file=$3			#output file basename
	
	file_format=${in_file#*.}
	line_count=`wc -l $in_file | awk '{print $1}'`		#number of lines in input file
	seq_count=$(( $line_count / $lpseq ))				#number of sequences in input file
	seqpfile=$(( $seq_count / 10 ))						#base number of sequences per output file
	sup_seq=$(( $seq_count % 10 ))						#number of output files with one supplementary sequence
	
	for ((i=0; 10-$i; i++)); do
		if [ $i -eq 0 ]; then
			start_seq=1
			if [ $sup_seq -eq 0 ]; then
				stop_seq=$seqpfile
			else
				stop_seq=$(( $seqpfile+1 ))
				sup_seq=$(( $sup_seq-1 ))
			fi
		else
			start_seq=$(( $stop_seq+1 ))
			if [ $sup_seq -eq 0 ]; then
				stop_seq=$(( $stop_seq+$seqpfile ))
			else
				stop_seq=$(( $stop_seq+$seqpfile+1 ))
				sup_seq=$(( $sup_seq-1 ))
			fi
		fi
		
		start_index=$(( $start_seq*$lpseq-$lpseq+1 ))
		stop_index=$(( $stop_seq*$lpseq ))
		
		sed -n ''$start_index','$stop_index'p' $in_file > "$out_file$i.$file_format"
	done
}
###############################################################################################################


#Filling of the training and testing datasets (Ti, Fi, Bti, Bfi + methylated, i from 0 to 9)

echo "Filling foreground and background files for $folder"
line_count=`wc -l $peaks_fg_bed | awk '{print $1}'`
echo "$line_count lines/sequences to process for $folder"
start_time=`date +%s`

### !!! Take care of the number of lines per sequence in the different files and chop them accordingly !!! #####
### (1 for BED, 2 or more for fasta) ########################################################################### 

chopFile $peaks_fg_bed 1 $foreground/test/bed/F
chopFile $peaks_fg_seq 2 $foreground/test/fasta/F
chopFile $peaks_bg_bed 1 $background/test/bed/Bf
chopFile $peaks_bg_seq 2 $background/test/fasta/Bf


#Concatenating test files to create (not) corresponding train files: test 1-9 -> train 0, test 0,2-9 -> train 1 etc.
for ((i=0; 10-$i; i++)); do
	ls $foreground/test/bed/F[0-9].bed | grep -v F$i.bed | xargs cat > $foreground/train/bed/T$i.bed
	ls $foreground/test/fasta/F[0-9].fa | grep -v F$i.fa | xargs cat > $foreground/train/fasta/T$i.fa
	ls $background/test/bed/Bf[0-9].bed | grep -v Bf$i.bed | xargs cat > $background/train/bed/Bt$i.bed
	ls $background/test/fasta/Bf[0-9].fa | grep -v Bf$i.fa | xargs cat > $background/train/fasta/Bt$i.fa
done

if [ $methylation_FIMO = true ]; then
	chopFile $peaks_fg_methyl_seq 2 $foreground/test/methyl_fasta/methyl_F
	chopFile $peaks_bg_methyl_seq 2 $background/test/methyl_fasta/methyl_Bf

	for ((i=0; 10-$i; i++)); do
		ls $foreground/test/methyl_fasta/methyl_F[0-9].fa | grep -v F$i.fa | xargs cat > $foreground/train/methyl_fasta/methyl_T$i.fa
		ls $background/test/methyl_fasta/methyl_Bf[0-9].fa | grep -v Bf$i.fa | xargs cat > $background/train/methyl_fasta/methyl_Bt$i.fa
	done
fi

end_time=`date +%s`
runtime=$((end_time-start_time))
echo "Runtime = $runtime s for $folder"
