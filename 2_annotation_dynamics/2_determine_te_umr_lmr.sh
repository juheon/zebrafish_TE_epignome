#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --array=1-11
#SBATCH --mem=2G
#SBATCH --job-name=te_umr_lmr_tissue

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID


# TE classes
bed_te=danRer10.TE_frag.bed.gz


# INPUT DATA
list=tissues-e.txt
tissue=$( cat $list | sed "${ID}q;d" )

dir_in=0_cres
bed_umr=${dir_in}/umr/fylab_WGBS_zt_${tissue}_UMR.bed.gz
bed_lmr=${dir_in}/lmr/fylab_WGBS_zt_${tissue}_LMR.bed.gz


# OUTPUT
dir_out=2_dynamic
out=${dir_out}/danRer10.TEs_${tissue}_UMR_LMR.bed.gz
#out=${dir_out}/te_dynamic_umr_lmr_across_10tissues.txt


# COMMANDS
intersectBed -sorted -wao -a $bed_te -b $bed_umr $bed_lmr |
    groupBy -g 1,2,3,8 -c 12 -o sum |
    groupBy -c 5 -o max -full |
    awk 'BEGIN{OFS="\t"} $4==1 {$4="UMR"} $4==2 {$4="LMR"} {print $1,$2,$3,$4}' | gzip -nc > $out

