#!/bin/bash
# Author: Hyung Joo Lee

### SOFTWARE
ml samtools/1.9

workdir=[]

cpus=12


### INPUT
list=${workdir}/samples.txt


for ID in {1..26}
do

    # INPUT
    sample=$( cat $list | sed "${ID}q;d" )
    dir_in=${workdir}/2_star_2pass
    bam_in=${dir_in}/fylab_RNA_${sample}.Aligned.out.bam

    # Output files
    bam_out=${dir_in}/fylab_RNA_${sample}.Aligned.sorted.out.bam

    # COMMANDS
    samtools sort -m 5G -o $bam_out -T /tmp/$sample -@ $cpus $bam_in
    samtools index -@ $cpus $bam_out

done



