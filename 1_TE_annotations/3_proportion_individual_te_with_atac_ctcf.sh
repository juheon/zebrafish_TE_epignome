#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=2,8
#SBATCH --job-name=propor_ind_te

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID


# TE 6 bed files
te=danRer10.TE_frag.bed.gz
dna=danRer10.DNA_frag.bed.gz
ltr=danRer10.LTR_frag.bed.gz
line=danRer10.LINE_frag.bed.gz
sine=danRer10.SINE_frag.bed.gz
rc=danRer10.RC_frag.bed.gz

te_bed=$(ls danRer10.{TE,DNA,LTR,LINE,SINE,RC}_frag.bed.gz )


# INPUT DATA
list=tissues-e.txt
tissue=$( cat $list | sed "${ID}q;d" )

dir_in=0_cres
tad=${dir_in}/TADs/TAD_boundaries_${tissue}.atac_ctcf.bed.gz
loop=${dir_in}/loop_anchors/loop_anchors_${tissue}.atac_ctcf.bed.gz


# OUTPUT
dir_out=2_dynamic
out=${dir_out}/fraction_atac_ctcf_of_TE_in_${tissue}.txt


# COMMANDS

for bed in $te_bed
do

	name=${bed##*/}
	name=${name%.bed.gz}
	name=${name#danRer10.}

	cnt_te_tad=$( intersectBed -u -a $bed -b $tad | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_tad=$( awk -v overlap="$cnt_te_tad" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_tad}\t${cnt_te_tad}\t${cnt_te}\tTAD_boundary\t${name}\t$tissue"

	cnt_te_loop=$( intersectBed -u -a $bed -b $loop | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_loop=$( awk -v overlap="$cnt_te_loop" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_loop}\t${cnt_te_loop}\t${cnt_te}\tLoop_anchor\t${name}\t$tissue"

done > $out
