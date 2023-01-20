#!/usr/bin/env bash

if [ "$1" = "--getorf_only" ] || [ "$1" = "-g" ]
then

	#ARGS
	reference_genome=$2
	min_orf_size=$3
	translation_table=$4
	force_svg=$5

	if [ -z "${reference_genome}" ] || [ -z "${min_orf_size}" ] || [ -z "${translation_table}" ] || [ -z "${force_svg}" ]
	then 
			echo """
codingDiv.sh v1.0 

-g, --getorf_only
    Produce only SVG plot of protein predictions, useful when no microdiversity data available

Positional arguments: 
1- Reference genome / Studied genome (FASTA)

2- Minimal ORF size (in nucleotides) [integer]

3- Translation table number used by EMBOSS getorf - https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi [integer 1-23]

4- Force SVG for a very large genome, over 100 kilobases [Y|N]

This last option is not recomended as it will generate a very large SVG file.
A better option would be splitting your genomes in several regions.

Cite us:

CodingDiv : visualize SNP-level microdiversity to discriminate between coding and noncoding regions.
Eric Olo Ndela & François Enault (2023, unpublished).
Laboratoire Microorganismes Genome & Environnement (LMGE)
Clermont-Auvergne University (UCA)
		        """
	else


		#Fetching reference name and size
		genome_name=$(grep '>' $reference_genome |awk '{print $1}' |sed -re 's/>//')
		genome_size=$(grep -v '>' $reference_genome | sed -z 's/\n//g' |wc -c)

		#file name without dir and extension
		file_name=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))

		#####################################################
		#FILE NAMES

		ref_orfs=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_orfs.faa

		ref_orfs_tsv=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_orfs.tsv

		ref_trslt=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_full.faa

		#Prodigal & phanotate
		prodigal_faa=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_prodigal.faa

		prodigal_gbk=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_prodigal.gbk

		phanotate_tsv=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_phanotate.tsv

		######################################################

		#Protein prediction
		echo "Protein prediction..."
		{
		printf "\n\n### getorf -sequence $reference_genome -outseq $ref_orfs -table $translation_table -minsize $min_orf_size -find 1 -circular N -reverse Y\n\n"
		getorf -sequence $reference_genome -outseq $ref_orfs -table $translation_table -minsize $min_orf_size -find 1 -circular N -reverse Y
		exit_code21=$?

		grep '>' $ref_orfs | awk '{if($5=="(REVERSE") {print $1"-\t"$2"\t"$4} else {print $1"+\t"$2"\t"$4}}' |sed -re 's/>// ; s/\[// ; s/\]//' > $ref_orfs_tsv

		####################prodigal & phanotate

		printf "\n\n### prodigal -i $reference_genome -o $prodigal_gbk -a $prodigal_faa -p meta\n\n"
		prodigal -i $reference_genome -o $prodigal_gbk -a $prodigal_faa -p meta
		exit_code22=$?

		printf "\n\n### phanotate.py $reference_genome > $phanotate_tsv\n\n"
		phanotate.py $reference_genome > $phanotate_tsv
		exit_code23=$?

		} &>>stdout.txt

		if [ $exit_code21 -eq 0 ] || [ $exit_code22 -eq 0 ] || [ $exit_code23 -eq 0 ]
		then

			if [ $genome_size -lt 100000 ]
			then

				echo "Plotting..."

				{
				########################## prodigal & phanotate maps
				printf "\n\n### prodigal_map.py $genome_name $genome_size $prodigal_faa\n\n"
				prodigal_map.py $genome_name $genome_size $prodigal_faa
				exit_code11=$?

				printf "\n\n### phanotate_map.py $genome_name $genome_size $phanotate_tsv\n\n"
				phanotate_map.py $genome_name $genome_size $phanotate_tsv
				exit_code12=$?

				printf "\n\n### getorf_map_pos.py $genome_name $genome_size $ref_orfs\n\n"
				getorf_map_pos.py $genome_name $genome_size $ref_orfs
				exit_code13=$?

				printf "\n\n### getorf_map_neg.py $genome_name $genome_size $ref_orfs\n\n"
				getorf_map_neg.py $genome_name $genome_size $ref_orfs
				exit_code14=$?

				printf "\n\n### svg_stack.py --direction=V --margin=15 $genome_name'_prodigal.svg' $genome_name'_phanotate.svg' $genome_name'_pnps.svg' pnps_legend.svg $genome_name'_neg_strand_pnps.svg'  $genome_name'_bar_chart.svg' > $file_name'_predictions.svg'\n\n"
				svg_stack.py --direction=V --margin=15 $genome_name"_prodigal.svg" $genome_name"_phanotate.svg" $genome_name"_getorf_pos.svg" $genome_name"_getorf_neg.svg" > $file_name"_predictions.svg"
				exit_code15=$?

				rm $genome_name"_prodigal.svg" $genome_name"_phanotate.svg" $genome_name"_getorf_pos.svg" $genome_name"_getorf_neg.svg"

				} &>>stdout.txt

				if [ $exit_code11 -eq 0 ] || [ $exit_code12 -eq 0 ] || [ $exit_code13 -eq 0 ] || [ $exit_code14 -eq 0 ] || [ $exit_code15 -eq 0 ]
				then	

					echo "DONE"
					echo "Genomic maps are plotted into "$file_name"_predictions.svg"
					echo "prodigal : " $prodigal_faa "and" $prodigal_gbk
					echo "phanotate :" $phanotate_tsv
					echo "getorf :" $ref_orfs

				else
					echo "#################################################"
					echo "There was an error, check stdout.txt for details"
					echo "#################################################"	
				fi

			elif [ $genome_size -gt 100000 ] && [ $force_svg = "Y" ]

				echo "Plotting..."

				{
				########################## prodigal & phanotate maps
				printf "\n\n### prodigal_map.py $genome_name $genome_size $prodigal_faa\n\n"
				prodigal_map.py $genome_name $genome_size $prodigal_faa
				exit_code11=$?

				printf "\n\n### phanotate_map.py $genome_name $genome_size $phanotate_tsv\n\n"
				phanotate_map.py $genome_name $genome_size $phanotate_tsv
				exit_code12=$?

				printf "\n\n### getorf_map_pos.py $genome_name $genome_size $ref_orfs\n\n"
				getorf_map_pos.py $genome_name $genome_size $ref_orfs
				exit_code13=$?

				printf "\n\n### getorf_map_neg.py $genome_name $genome_size $ref_orfs\n\n"
				getorf_map_neg.py $genome_name $genome_size $ref_orfs
				exit_code14=$?

				printf "\n\n### svg_stack.py --direction=V --margin=15 $genome_name'_prodigal.svg' $genome_name'_phanotate.svg' $genome_name'_pnps.svg' pnps_legend.svg $genome_name'_neg_strand_pnps.svg'  $genome_name'_bar_chart.svg' > $file_name'_predictions.svg'\n\n"
				svg_stack.py --direction=V --margin=15 $genome_name"_prodigal.svg" $genome_name"_phanotate.svg" $genome_name"_getorf_pos.svg" $genome_name"_getorf_neg.svg" > $file_name"_predictions.svg"
				exit_code15=$?

				rm $genome_name"_prodigal.svg" $genome_name"_phanotate.svg" $genome_name"_getorf_pos.svg" $genome_name"_getorf_neg.svg"

				} &>>stdout.txt

				if [ $exit_code11 -eq 0 ] || [ $exit_code12 -eq 0 ] || [ $exit_code13 -eq 0 ] || [ $exit_code14 -eq 0 ] || [ $exit_code15 -eq 0 ]
				then	

					echo "DONE"
					echo "Genomic maps are plotted into "$file_name"_predictions.svg"
					echo "prodigal : " $prodigal_faa "and" $prodigal_gbk
					echo "phanotate :" $phanotate_tsv
					echo "getorf :" $ref_orfs

				else
					echo "#################################################"
					echo "There was an error, check stdout.txt for details"
					echo "#################################################"	
				fi

			else
				echo "No SVG was produced due to the large size of the genome"
				echo "DONE"
				echo "prodigal : " $prodigal_faa "and" $prodigal_gbk
				echo "phanotate :" $phanotate_tsv
				echo "getorf :" $ref_orfs
			fi

		else
			echo "#################################################"
			echo "There was an error, check stdout.txt for details"
			echo "#################################################"
		fi
	fi

else

	#ARGS
	reference_genome=$1

	reads_contigs=$2

	#To manage globs
	if [[ $reads_contigs == *"*"* ]]; then
	  reads_contigs=$(awk 'gsub(/\\/,"",$0)' <(echo "$reads_contigs"))
	fi


	min_orf_size=$3

	translation_table=$4

	min_reads=$5

	min_percentage=$6

	#for mapping steps
	num_threads=$7

	force_svg="$8"

	#Fetching reference name and size
	genome_name=$(grep '>' $reference_genome |awk '{print $1}' |sed -re 's/>//')
	genome_size=$(grep -v '>' $reference_genome | sed -z 's/\n//g' |wc -c)

	#deletng stdout.txt file from previous run 
	{
	rm stdout.txt
	} &>>stdout.txt

	if [ -z "${reference_genome}" ] || [ -z "${reads_contigs}" ] || [ -z "${min_orf_size}" ] || [ -z "${translation_table}" ] || [ -z "${min_reads}" ] || [ -z "${min_percentage}" ] || [ -z "${num_threads}" ]
	then 
		echo """
	codingDiv.sh v1.0 

	Tags:
	-g, --getorf_only
	    Produce only SVG plot of protein predictions, useful when no microdiversity data available

	Positional arguments: 
	1- Reference genome / Studied genome (FASTA)

	2- Reads or contigs to map (FASTA or FASTQ)
	Here you can input a set of metagenomes with a glob, escaping the character as follows: \"\*.fastq.gz\"

	3- Minimal ORF size (in nucleotides) [integer]

	4- Translation table number used by EMBOSS getorf - https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi [integer 1-23]

	5- Minimal number of reads required to evaluate an SNP [integer]

	6- Minimal % of total mapped reads required to evaluate an SNP [double]

	7- Number of CPU allowed for mapping [integer]

	8- Force SVG for a very large genome, over 100 kilobases [Y|N]

	This last option is not recomended as it will generate a very large SVG file.
	A better option would be splitting your genomes in several regions.

	Cite us:

	CodingDiv : visualize SNP-level microdiversity to discriminate between coding and noncoding regions.
	Eric Olo Ndela & François Enault (2023, unpublished).
	Laboratoire Microorganismes Genome & Environnement (LMGE)
	Clermont-Auvergne University (UCA)
	        """
	else
		#file name without dir and extension
		file_name=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))

		if [ ! -f "$reference_genome" ]
		then 
			echo "ERROR : One of the FASTA/FASTQ files does not exist"
		else 

			######################################################################
			#FILE NAMES

			ref_orfs=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_orfs.faa

			ref_orfs_tsv=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_orfs.tsv

			ref_trslt=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_full.faa

			snp_file=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome))).final_variants_with_depth.tsv

			codon_file=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_codon_freq.tsv

			depth_file=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome))).depth

			#########################################################################
			#Prodigal & phanotate
			prodigal_faa=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_prodigal.faa

			prodigal_gbk=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_prodigal.gbk

			phanotate_tsv=$(awk -F'.' '{print $1}' <(echo $(basename $reference_genome)))_phanotate.tsv

			#########################################################################
			echo "Protein prediction..."

			{

			mkdir mapping_files prot_prediction output_files

			####################orf prediction
			printf "\n\n### getorf -sequence $reference_genome -outseq $ref_orfs -table $translation_table -minsize $min_orf_size -find 1 -circular N -reverse Y\n\n"
			getorf -sequence $reference_genome -outseq $ref_orfs -table $translation_table -minsize $min_orf_size -find 1 -circular N -reverse Y

			grep '>' $ref_orfs | awk '{if($5=="(REVERSE") {print $1"-\t"$2"\t"$4} else {print $1"+\t"$2"\t"$4}}' |sed -re 's/>// ; s/\[// ; s/\]//' > $ref_orfs_tsv

			####################prodigal & phanotate

			printf "\n\n### prodigal -i $reference_genome -o $prodigal_gbk -a $prodigal_faa -p meta\n\n"
			prodigal -i $reference_genome -o $prodigal_gbk -a $prodigal_faa -p meta

			printf "\n\n### phanotate.py $reference_genome > $phanotate_tsv\n\n"
			phanotate.py $reference_genome > $phanotate_tsv

			#####################prediction full reading frames
			printf "\n\n### translate.py $reference_genome > $ref_trslt\n\n"
			translate.py $reference_genome > $ref_trslt

			} &>>stdout.txt

			#####################Mapping
			echo "Mapping... "$reads_contigs" to "$1

			{
			printf "\n\n### map_var_call.sh $reference_genome $2 $num_threads\n\n"
			map_var_call.sh $reference_genome $2 $num_threads
			exit_code=$?
			} &>>stdout.txt

			if [ $exit_code -eq 0 ] && [[ $(awk '{a+=length($NF)} END {print a}' $snp_file) -ne $genome_size ]] && [ -s $depth_file ]
			then

				########################SNP assessment, on full length reading frames
				echo "Checking effect of SNPs on protein sequence..."

				{
				printf "\n\n### snp_check_full_length_ORFs.py $ref_trslt $snp_file $reference_genome $translation_table > $ref_trslt.snp_summary.tsv\n\n"
				snp_check_full_length_ORFs.py $ref_trslt $snp_file $reference_genome $translation_table > $ref_trslt.snp_summary.tsv
				exit_code1=$?

				#on ORFs from start to stop (both strands)
				printf "\n\n### snp_check_all_prot_REVERSE.py $ref_orfs $snp_file $reference_genome $translation_table > $ref_orfs.snp_summary.tsv\n\n"
				snp_check_all_prot_REVERSE.py $ref_orfs $snp_file $reference_genome $translation_table > $ref_orfs.snp_summary.tsv
				exit_code2=$?
				
				} &>>stdout.txt

				if [ $exit_code1 -eq 0 ] && [ $exit_code2 -eq 0 ] 
				then

					echo "Producing plots & moving files to their directories..."

					{


					############################ORFs pNpS
					#name of the genome and size /. full orfs / getorf
					printf "\n\n### pnps.R $genome_name $genome_size $ref_orfs.snp_summary.tsv $min_reads $min_percentage $ref_orfs_tsv\n\n"
					pnps.R $genome_name $genome_size $ref_orfs.snp_summary.tsv $min_reads $min_percentage $ref_orfs_tsv
					exit_code1=$?

					#Producing plots
					if [ $genome_size -lt 100000 ]
					then

						##########################MAPPING plots (full orfs of course)
						printf "\n\n### snp_plot.R $genome_name $genome_size $ref_trslt.snp_summary.tsv $min_reads $min_percentage\n\n"
						snp_plot.R $genome_name $genome_size $ref_trslt.snp_summary.tsv $min_reads $min_percentage
						exit_code2=$?

						printf "\n\n### map_with_pnps.py $genome_name $genome_size $ref_orfs\n\n"
						map_with_pnps.py $genome_name $genome_size $ref_orfs
						exit_code3=$?

						printf "\n\n### map_with_pnps_neg_strand.py $genome_name $genome_size $ref_orfs $depth_file\n\n"
						map_with_pnps_neg_strand.py $genome_name $genome_size $ref_orfs $depth_file
						exit_code4=$?

						########################## prodigal & phanotate maps
						printf "\n\n### prodigal_map.py $genome_name $genome_size $prodigal_faa\n\n"
						prodigal_map.py $genome_name $genome_size $prodigal_faa
						exit_code5=$?

						printf "\n\n### phanotate_map.py $genome_name $genome_size $phanotate_tsv\n\n"
						phanotate_map.py $genome_name $genome_size $phanotate_tsv
						exit_code6=$?

					elif [ $genome_size -gt 100000 ] && [ $force_svg = "Y" ]
					then
						##########################MAPPING plots (full orfs of course)
						printf "\n\n### snp_plot.R $genome_name $genome_size $ref_trslt.snp_summary.tsv $min_reads $min_percentage\n\n"
						snp_plot.R $genome_name $genome_size $ref_trslt.snp_summary.tsv $min_reads $min_percentage
						exit_code2=$?

						printf "\n\n### map_with_pnps.py $genome_name $genome_size $ref_orfs\n\n"
						map_with_pnps.py $genome_name $genome_size $ref_orfs
						exit_code3=$?

						printf "\n\n### map_with_pnps_neg_strand.py $genome_name $genome_size $ref_orfs $depth_file\n\n"
						map_with_pnps_neg_strand.py $genome_name $genome_size $ref_orfs $depth_file
						exit_code4=$?

						########################## prodigal & phanotate maps
						printf "\n\n### prodigal_map.py $genome_name $genome_size $prodigal_faa\n\n"
						prodigal_map.py $genome_name $genome_size $prodigal_faa
						exit_code5=$?

						printf "\n\n### phanotate_map.py $genome_name $genome_size $phanotate_tsv\n\n"
						phanotate_map.py $genome_name $genome_size $phanotate_tsv
						exit_code6=$?

					fi

					} &>>stdout.txt

					if [ $exit_code1 -eq 0 ] && [ $exit_code2 -eq 0 ] && [ $exit_code3 -eq 0 ] && [ $exit_code4 -eq 0 ] && [ $exit_code5 -eq 0 ] && [ $exit_code6 -eq 0 ]
					then

						{
						########################## moving files 

						mv $genome_name"_pnps.svg" output_files/

						mv $genome_name"_neg_strand_pnps.svg" output_files/

						mv $genome_name"_prodigal.svg" output_files/

						mv $genome_name"_phanotate.svg" output_files/

						mv $reference_genome.* mapping_files/

						mv pnps_legend.svg output_files/

						mv *.faa prot_prediction/

						mv *.gbk prot_prediction/

						mv $phanotate_tsv prot_prediction/

						mv $snp_file output_files/


						mv *.bam.* mapping_files/

						mv *.bam mapping_files/

						mv *.var.* mapping_files/

						mv *.bcf mapping_files/

						mv *.depth mapping_files/

						mv Rplots.pdf mapping_files/


						mv $ref_trslt.snp_summary.tsv output_files/

						mv $ref_orfs.snp_summary.tsv output_files/

						mv $ref_orfs_tsv prot_prediction/

						mv $genome_name"_pnps.tsv" output_files/

						mv $genome_name"_bar_chart.svg" output_files/

						cd output_files

						if [ $genome_size -lt 100000 ]
						then
							printf "\n\n### svg_stack.py --direction=V --margin=15 $genome_name'_prodigal.svg' $genome_name'_phanotate.svg' $genome_name'_pnps.svg' pnps_legend.svg $genome_name'_neg_strand_pnps.svg'  $genome_name'_bar_chart.svg' > $file_name'_summary.svg'\n\n"
							svg_stack.py --direction=V --margin=15 $genome_name"_prodigal.svg" $genome_name"_phanotate.svg" $genome_name"_pnps.svg" pnps_legend.svg $genome_name"_neg_strand_pnps.svg"  $genome_name"_bar_chart.svg" > $file_name"_summary.svg"
							exit_code=$?

						elif [ $genome_size -gt 100000 ] && [ $force_svg = "Y" ]
						then
							printf "\n\n### svg_stack.py --direction=V --margin=15 $genome_name'_prodigal.svg' $genome_name'_phanotate.svg' $genome_name'_pnps.svg' pnps_legend.svg $genome_name'_neg_strand_pnps.svg'  $genome_name'_bar_chart.svg' > $file_name'_summary.svg'\n\n"
							svg_stack.py --direction=V --margin=15 $genome_name"_prodigal.svg" $genome_name"_phanotate.svg" $genome_name"_pnps.svg" pnps_legend.svg $genome_name"_neg_strand_pnps.svg"  $genome_name"_bar_chart.svg" > $file_name"_summary.svg"
							exit_code=$?
						fi

						#Avoiding black background

						sed -i 's/#333333;"/#333333; fill: none;"/g' $file_name"_summary.svg"

						#avoiding black squares around maps

						sed -i 's/"fill:#ffffff;"\/>/"fill:#ffffff;stroke:#ffffff"\/>/g' $file_name"_summary.svg"

						#avoiding white text

						sed -i 's/fill: none;$//g' $file_name"_summary.svg"

						mkdir ../final_results

						mv ../summary_table.tsv ../final_results/$file_name"_summary_table.tsv"

						mv $file_name"_summary.svg" ../final_results/

						mv ../full_snp_table.tsv ../final_results/$file_name"_full_snp_table.tsv"

						} &>>stdout.txt

						echo "-- ALL DONE check out the final_results/ dir --"

					else

						{
						########################## moving files 

						mv $genome_name"_pnps.svg" output_files/

						mv $genome_name"_neg_strand_pnps.svg" output_files/

						mv $genome_name"_prodigal.svg" output_files/

						mv $genome_name"_phanotate.svg" output_files/

						mv $reference_genome.* mapping_files/

						mv pnps_legend.svg output_files/

						mv *.faa prot_prediction/

						mv *.gbk prot_prediction/

						mv $phanotate_tsv prot_prediction/

						mv $snp_file output_files/


						mv *.bam.* mapping_files/

						mv *.bam mapping_files/

						mv *.var.* mapping_files/

						mv *.bcf mapping_files/

						mv *.depth mapping_files/

						mv Rplots.pdf mapping_files/


						mv $ref_trslt.snp_summary.tsv output_files/

						mv $ref_orfs.snp_summary.tsv output_files/

						mv $ref_orfs_tsv prot_prediction/

						mv $genome_name"_pnps.tsv" output_files/

						mv $genome_name"_bar_chart.svg" output_files/

						cd output_files

						if [ $genome_size -lt 100000 ]
						then
							printf "\n\n### svg_stack.py --direction=V --margin=15 $genome_name'_prodigal.svg' $genome_name'_phanotate.svg' $genome_name'_pnps.svg' pnps_legend.svg $genome_name'_neg_strand_pnps.svg'  $genome_name'_bar_chart.svg' > $file_name'_summary.svg'\n\n"
							svg_stack.py --direction=V --margin=15 $genome_name"_prodigal.svg" $genome_name"_phanotate.svg" $genome_name"_pnps.svg" pnps_legend.svg $genome_name"_neg_strand_pnps.svg"  $genome_name"_bar_chart.svg" > $file_name"_summary.svg"
							exit_code=$?

						elif [ $genome_size -gt 100000 ] && [ $force_svg = "Y" ]
						then
							printf "\n\n### svg_stack.py --direction=V --margin=15 $genome_name'_prodigal.svg' $genome_name'_phanotate.svg' $genome_name'_pnps.svg' pnps_legend.svg $genome_name'_neg_strand_pnps.svg'  $genome_name'_bar_chart.svg' > $file_name'_summary.svg'\n\n"
							svg_stack.py --direction=V --margin=15 $genome_name"_prodigal.svg" $genome_name"_phanotate.svg" $genome_name"_pnps.svg" pnps_legend.svg $genome_name"_neg_strand_pnps.svg"  $genome_name"_bar_chart.svg" > $file_name"_summary.svg"
							exit_code=$?
						fi

						#Avoiding black background

						sed -i 's/#333333;"/#333333; fill: none;"/g' $file_name"_summary.svg"

						#avoiding black squares around maps

						sed -i 's/"fill:#ffffff;"\/>/"fill:#ffffff;stroke:#ffffff"\/>/g' $file_name"_summary.svg"

						#avoiding white text

						sed -i 's/fill: none;$//g' $file_name"_summary.svg"

						mkdir ../final_results

						mv ../summary_table.tsv ../final_results/$file_name"_summary_table.tsv"

						mv $file_name"_summary.svg" ../final_results/

						mv ../full_snp_table.tsv ../final_results/$file_name"_full_snp_table.tsv"

						} &>>stdout.txt
			 
						echo "#################################################"
						echo "There was an error, check stdout.txt for details"
						echo "#################################################"
					fi

				else
					echo "#################################################"
					echo "There was an error, check stdout.txt for details"
					echo "#################################################"
				fi

			elif [ $exit_code -eq 0 ] && [[ $(awk '{a+=length($NF)} END {print a}' $snp_file) == $genome_size ]] || [ ! -s $depth_file ]
			then
				echo "#####################################"
				echo "No SNPs were found, execution halted"
				echo "#####################################"
			else
				echo "#################################################"
				echo "There was an error, check stdout.txt for details"
				echo "#################################################"
			fi
		fi
	fi

	{
	########################## moving files even if errors

	mv $genome_name"_pnps.svg" output_files/

	mv $genome_name"_neg_strand_pnps.svg" output_files/

	mv $genome_name"_prodigal.svg" output_files/

	mv $genome_name"_phanotate.svg" output_files/

	mv $reference_genome.* mapping_files/

	mv pnps_legend.svg output_files/

	mv *.faa prot_prediction/

	mv *.gbk prot_prediction/

	mv $phanotate_tsv prot_prediction/

	mv $snp_file output_files/


	mv *.bam.* mapping_files/

	mv *.bam mapping_files/

	mv *.var.* mapping_files/

	mv *.bcf mapping_files/

	mv *.depth mapping_files/

	mv Rplots.pdf mapping_files/


	mv $ref_trslt.snp_summary.tsv output_files/

	mv $ref_orfs.snp_summary.tsv output_files/

	mv $ref_orfs_tsv prot_prediction/

	mv $genome_name"_pnps.tsv" output_files/

	mv $genome_name"_bar_chart.svg" output_files/

	cd output_files

	#Avoiding black background

	sed -i 's/#333333;"/#333333; fill: none;"/g' $file_name"_summary.svg"

	#avoiding black squares around maps

	sed -i 's/"fill:#ffffff;"\/>/"fill:#ffffff;stroke:#ffffff"\/>/g' $file_name"_summary.svg"

	#avoiding white text

	sed -i 's/fill: none;$//g' $file_name"_summary.svg"

	mkdir ../final_results

	mv ../summary_table.tsv ../final_results/$file_name"_summary_table.tsv"

	mv $file_name"_summary.svg" ../final_results/

	mv ../full_snp_table.tsv ../final_results/$file_name"_full_snp_table.tsv"

	} &>>stdout.txt
fi
