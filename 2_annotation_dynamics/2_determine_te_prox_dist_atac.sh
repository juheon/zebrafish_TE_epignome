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
bed_prox=${dir_in}/proximalatacseq/${tissue}.proximal.open.bed
bed_dist=${dir_in}/distalatacseq/${tissue}.distal.open.bed


# OUTPUT
dir_out=2_dynamic
out=${dir_out}/danRer10.TEs_${tissue}_ATAC.bed.gz


# COMMANDS
intersectBed -sorted -wao -a $bed_te -b <(sort -k1,1 -k2,2n $bed_prox) <(sort -k1,1 -k2,2n $bed_dist) |
    groupBy -g 1,2,3,8 -c 13 -o sum |
    groupBy -c 5 -o max -full |
    awk 'BEGIN{OFS="\t"} $4==1 {$4="Proximal"} $4==2 {$4="Distal"} {print $1,$2,$3,$4}' | gzip -nc > $out

