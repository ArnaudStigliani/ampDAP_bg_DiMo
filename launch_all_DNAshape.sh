#!/usr/bin/env bash


>>../results/stderr.txt
>>../results/stdout.txt

path_scripts=/storage/mathelierarea/processed/arnaudst/ampDAP_bg_DiMO/scripts
data=/storage/mathelierarea/processed/arnaudst/data
results_shared=/storage/mathelierarea/processed/arnaudst/ampDAP_bg_DiMO/results
path_BiasAway=/lsc/BiasAway

genome=$data/Arabidopsis_thaliana.TAIR10.dna.genome.fa
methyl_genome=$data/AraTha_TAIR10_methyl_genome.fa
methyl_alphabet=$data/alphabet
chr_lengths=$data/AraTha_chr_length.genome


#ls -d $data/dap_data_v4/peaks/*/* | grep "colamp" | sed -e 's/_colamp.*//g' | xargs -I{} --max-proc=8 bash -c \


# ls -d $data/dap_data_v4/peaks/[N-Z]*/* | grep "colamp" | sed -e 's/_colamp.*//g' | xargs -I{} --max-proc=8 bash -c \
ls -d $data/dap_data_v4/peaks/C2H2_tnt/At5g66730_colamp_a | grep "colamp" | sed -e 's/_colamp.*//g' | xargs -I{} --max-proc=5 bash -c \
'
	path_scripts=/storage/mathelierarea/processed/arnaudst/ampDAP_bg_DiMO/scripts
	data=/storage/mathelierarea/processed/arnaudst/data
	results_shared=/storage/mathelierarea/processed/arnaudst/ampDAP_bg_DiMO/results
	results_DAP=$results_shared/DAP_DNAshape

	genome=$data/Arabidopsis_thaliana.TAIR10.dna.genome.fa
	methyl_genome=$data/AraTha_TAIR10_methyl_genome.fa
	methyl_alphabet=$data/alphabet
	chr_lengths=$data/AraTha_chr_length.genome
	
	# genome=$0
	# methyl_genome=$1
	# methyl_alphabet=$2
	# chr_lengths=$3
	
	fg_folder=`ls -d $data/dap_data_v4/peaks/*/* | grep {} | grep -v "colamp"`
	bg_folder=`ls -d $data/dap_data_v4/peaks/*/* | grep {} | grep "colamp"`
	
	TF_family=$(basename $(dirname {}))
	TF_name=$(basename {})


	if [ -z "$fg_folder" ]; then 
		echo "The DAP-seq peaks are not available for $TF_name from $TF_family"; exit 1
	elif [ -z "$bg_folder" ]; then 
		echo "The ampDAP-seq peaks are not available for $TF_name from $TF_family"; exit 1
	fi

	peaks_fg=$fg_folder/"chr1-5"/"chr1-5_GEM_events.narrowPeak"
	peaks_bg=$bg_folder/"chr1-5"/"chr1-5_GEM_events.narrowPeak"
	
	nb_peaks=`wc -l $peaks_fg | cut -d " " -f1`
		
	if [ $nb_peaks -ge 1800 ]; then
		echo "$nb_peaks peaks for $TF_name from $TF_family, starting to compute"
		
		bash $path_scripts/DNAshapedTFBS_DAPseq_FIMO.sh -f $TF_family -n $TF_name -p $peaks_fg -b $peaks_bg -g $genome -c $methyl_genome -a $methyl_alphabet -y -q -t -l $chr_lengths -r $results_DAP
	
	# 	echo "TF $TF_name from $TF_family w/ methylation for hits and shapes from methylated motif: Done"
	
	# 	bash $path_scripts/DNAshapedTFBS_DAPseq_FIMO.sh -f $TF_family -n $TF_name -p $peaks_fg -b $peaks_bg -g $genome -c $methyl_genome -a $methyl_alphabet -y -q -l $chr_lengths -r $results_DAP
	
	# 	echo "TF $TF_name from $TF_family w/ methylation for hits only (not for shapes) from methylated motif: Done"
	
	# 	bash $path_scripts/DNAshapedTFBS_DAPseq_FIMO.sh -f $TF_family -n $TF_name -p $peaks_fg -b $peaks_bg -g $genome -c $methyl_genome -a $methyl_alphabet -y -t -l $chr_lengths -r $results_DAP
	
	# 	echo "TF $TF_name from $TF_family w/ methylation for shapes only (not for hits) from methylated motif: Done"
	
	# 	bash $path_scripts/DNAshapedTFBS_DAPseq_FIMO.sh -f $TF_family -n $TF_name -p $peaks_fg -b $peaks_bg -g $genome -c $methyl_genome -a $methyl_alphabet -y -l $chr_lengths -r $results_DAP
	
	# 	echo "TF $TF_name from $TF_family w/o methylation from methylated motif: Done"
		
	# 	bash $path_scripts/DNAshapedTFBS_DAPseq_FIMO.sh -f $TF_family -n $TF_name -p $peaks_fg -b $peaks_bg -g $genome -c $methyl_genome -a $methyl_alphabet -q -t -l $chr_lengths -r $results_DAP
	
	# 	echo "TF $TF_name from $TF_family w/ methylation for hits and shapes from non methylated motif: Done"
	
	# 	bash $path_scripts/DNAshapedTFBS_DAPseq_FIMO.sh -f $TF_family -n $TF_name -p $peaks_fg -b $peaks_bg -g $genome -c $methyl_genome -a $methyl_alphabet -q -l $chr_lengths -r $results_DAP
	
	# 	echo "TF $TF_name from $TF_family w/ methylation for hits only (not for shapes) from non methylated motif: Done"
	
	# 	bash $path_scripts/DNAshapedTFBS_DAPseq_FIMO.sh -f $TF_family -n $TF_name -p $peaks_fg -b $peaks_bg -g $genome -c $methyl_genome -a $methyl_alphabet -t -l $chr_lengths -r $results_DAP
	
	# 	echo "TF $TF_name from $TF_family w/ methylation for shapes only (not for hits) from non methylated motif: Done"
	
	# 	bash $path_scripts/DNAshapedTFBS_DAPseq_FIMO.sh -f $TF_family -n $TF_name -p $peaks_fg -b $peaks_bg -g $genome -l $chr_lengths -r $results_DAP
	
	# 	echo "TF $TF_name from $TF_family w/o methylation from non methylated motif: Done"
	
	# 	python $path_scripts/plot_hits_distribution.py $results_DAP/$TF_family/$TF_name
	
	else 
	        echo "$nb_peaks peaks for $TF_name from $TF_family, not enough to compute"
	fi
' 
#  $genome $methyl_genome $methyl_alphabet $chr_lengths -- {} 1>>$results_shared/stdout.txt 2>>$results_shared/stderr.txt

#python $path_scripts/plot_graphs.py $results_shared/DAP_DNAshape/AUPRC_scores.txt $results_shared/DAP_DNAshape/AUROC_scores.txt $data/methyl_sensitivity.tsv
