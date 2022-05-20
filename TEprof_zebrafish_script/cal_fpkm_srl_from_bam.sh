#!/bin/bash
# Author: Hyung Joo Lee

# SOFTWARE
ml samtools/1.9
ml bedtools/2.27.1

# INPUT and parameters
workdir=/scratch/hlee/fylab/rna/teprof2

list=${workdir}/samples.txt

gene="srl"
gene_bed_c=/scratch/hlee/SRA/${gene}_canonical.bed
gene_bed_te=/scratch/hlee/SRA/${gene}_te_derived.bed
coord="chr3:55950000-55970000"


# OUTPUT
out=${workdir}/count_bam_${gene}.txt


# COMMANDS
l_c=$( awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' $gene_bed_c )
l_te=$(awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' $gene_bed_te )

#for id in {2..90}
for id in {1..26}
do
    sample=$( cat $list | sed "${id}q;d" )
    input=${workdir}/2_star_2pass/${sample}.Aligned.sorted.out.bam
    log_star=${workdir}/2_star_2pass/${sample}.Log.final.out
    

    c_total=$( cat $log_star | grep "Uniquely mapped reads number" | cut -f2 )
    c_c=$( samtools view -b $input $coord |
        coverageBed -s -counts -sorted -split -a $gene_bed_c -b stdin |
        awk 'BEGIN{s=0} {s+=$7} END{print s}' )
    c_te=$( samtools view -b $input $coord |
        coverageBed -s -counts -sorted -split -a $gene_bed_te -b stdin |
        awk 'BEGIN{s=0} {s+=$7} END{print s}' )

    echo -e "$gene\t$sample\t$c_total\t$c_c\t$c_te\t$l_c\t$l_te" |
        awk -F"\t" -vOFS="\t" '{print $0,sprintf("%.4f\t%.4f",($4/($6/1000))/($3/1000000), ($5/($7/1000))/($3/1000000)) }'

done >> $out

