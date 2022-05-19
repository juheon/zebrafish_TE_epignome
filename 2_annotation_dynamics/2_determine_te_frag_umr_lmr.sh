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
tmp=${dir_out}/danRer10.TEs_${tissue}_UMR_LMR.tmp.bed
out=${dir_out}/danRer10.TEs_${tissue}_UMR_LMR.bed
#out=${dir_out}/te_dynamic_umr_lmr_across_10tissues.txt


# COMMANDS
intersectBed -sorted -wao -a $bed_te -b $bed_umr $bed_lmr |
    groupBy -g 1,2,3,4,5,6,7,8 -c 12 -o sum |
    groupBy -g 1,2,3,4,5,6,7 -c 9 -o max > $tmp

intersectBed -sorted -wao -a $bed_te -b $bed_umr $bed_lmr |
    groupBy -g 1,2,3,4,5,6,7,8 -c 12 -o sum |
    intersectBed -sorted -wao -f 1 -F 1 -a stdin -b $tmp | awk '$9==$(NF-1)' |
    groupBy -g 1,2,3,4,5,6,7 -c 8 -o first |
    awk 'BEGIN{OFS="\t"} $8==1 {$8="UMR"} $8==2 {$8="LMR"} {print $1,$2,$3,$4,$5,$6,$7,$8}' > $out

