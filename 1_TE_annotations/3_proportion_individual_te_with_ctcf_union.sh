#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
##SBATCH --array=1-11
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=propor_ind_te

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID

workdir=/scratch/twlab/hlee/zf_te


# TE bed file
dir_te=/scratch/twlab/hlee/genomes/danRer10/rmsk
te=${dir_te}/danRer10.TE.bed.gz
dna=${dir_te}/danRer10.DNA.bed.gz
ltr=${dir_te}/danRer10.LTR.bed.gz
line=${dir_te}/danRer10.LINE.bed.gz
sine=${dir_te}/danRer10.SINE.bed.gz
rc=${dir_te}/danRer10.RC.bed.gz

te_bed=$(ls ${dir_te}/danRer10.{TE,DNA,LTR,LINE,SINE,RC}_frag.bed.gz )


# INPUT DATA
list=${workdir}/tissues-e.txt
tissue=$( cat $list | sed "${ID}q;d" )

dir_in=/scratch/twlab/hlee/zf_te/0_cres/union
tad=${dir_in}/TAD_boundaries_union.ctcf.bed.gz
loop=${dir_in}/loop_anchors_union.ctcf.bed.gz


# OUTPUT
dir_out=${workdir}/2_dynamic
out=${dir_out}/fraction_ctcf_union_of_TE.txt


# COMMANDS

for bed in $te_bed
do

	name=${bed##*/}
	name=${name%.bed.gz}
	name=${name#danRer10.}

	cnt_te_tad=$( intersectBed -u -a $bed -b $tad | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_tad=$( awk -v overlap="$cnt_te_tad" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_tad}\t${cnt_te_tad}\t${cnt_te}\tTAD_boundary\t${name}\tUnion"

	cnt_te_loop=$( intersectBed -u -a $bed -b $loop | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_loop=$( awk -v overlap="$cnt_te_loop" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_loop}\t${cnt_te_loop}\t${cnt_te}\tLoop_anchor\t${name}\tUnion"

done > $out