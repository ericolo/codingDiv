#!/usr/bin/env bash


#MAPPING AND VARIANT CALLING

reference=$1

refname=$(awk -F'.' '{print $1}' <(echo $(basename $reference)))

reads=$2

#To manage globs
if [[ $reads == *"*"* ]]; then
  reads="${reads:1}"
fi

threads=$3

echo "Reference genome : $refname"
echo "Reads file/glob: $reads"

#MAPPING
bwa index $reference

bwa mem -t $threads $reference <(cat $reads) | samtools sort --threads $threads -o $refname.sorted.bam -

samtools index -@ $threads $refname.sorted.bam

echo "Mapping done. To view it, use the following command"
echo samtools tview $refname.sorted.bam $reference

#Checking the depth to adjust to variant calling 

samtools depth $refname.sorted.bam -o $refname.depth

depth=$(awk '{print $NF}' $refname.depth  |sort -rh |head -1)

#VARIANT CALLING
#We decided to skip indels because they do not tell us nothing in terms of codon constraints


bcftools mpileup -a FORMAT/AD --threads $threads --skip-indels --max-depth $depth -f $reference $refname.sorted.bam |bcftools call --threads $threads -m --skip-variants indels -p 3000 -A > $refname.variants.bcf



#In this last file DP = depth on the given position, and DP4= readsF,readsR,alt_readsF,alt_readsR, alt_reads will be 0,0 if no variant was found
#The AD tag (last column) is read count for each allele, in the order of the 5th column, comma separated 

grep -v '^##' $refname.variants.bcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$10}' |awk -F':' '{print $1"\t"$NF}' |sed -re 's/\t[0-9]+\/[0-9]+//' |grep -v '^#' > $refname.final_variants_with_depth.tsv


#step4 : selecting the high quality/frequent SNPs and creating files for the plots

awk '{if($6>=30) {print $2"\t"$5} else {print $2"\t""."}}' $refname.final_variants_with_depth.tsv |sed -re 's/\,//g' |awk '{if($2~/\./) {print $1"\t"1} else {print $1"\t"length($2)+1}}' > $refname.var.nuccount.tsv


awk '{if($6>=30) {print $2"\t"$NF} else {print $2"\t""0"}}' $refname.final_variants_with_depth.tsv  |sed -re 's/\,/\t/g' |awk '{if (NF==2) {print $1"\t""0"} else {print $1"\t"(($3+$4+$5)/($2+$3+$4+$5))*100} }' > $refname.var.read_percent.tsv


echo 'END'
echo "Variants were summarised in $refname.final_variants_with_depth.tsv"
echo "For more details check the raw file: $refname.variants.bcf"
echo "The last two columns represent total raw depth, and the number of reads accounting for variants, RESPCETIVELY"

echo ""

echo "High quality SNPs are written in plotting files : $refname.var.nuccount.tsv and $refname.var.read_percent.tsv"

echo ""

echo "If several nucs are possible at one column their depth is also printed in the last column of $refname.final_variants_with_depth.tsv"
