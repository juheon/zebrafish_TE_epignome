#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --array=1-11
#SBATCH --mem=2G
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=te_umr_lmr_tissue

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID


# TE classes
dir_te=/scratch/twlab/hlee/genomes/danRer10/rmsk
bed_te=${dir_te}/danRer10.TE_frag.bed.gz


# INPUT DATA
list=tissues-e.txt
tissue=$( cat $list | sed "${ID}q;d" )

dir_in=0_cres
bed_prox=${dir_in}/proximalatacseq/${tissue}.proximal.open.bed
bed_dist=${dir_in}/distalatacseq/${tissue}.distal.open.bed


# OUTPUT
dir_out=2_dynamic
tmp=${dir_out}/danRer10.TEs_${tissue}_ATAC.tmp.bed
out=${dir_out}/danRer10.TEs_${tissue}_ATAC.bed


# COMMANDS
intersectBed -sorted -wao -a $bed_te -b <(sort -k1,1 -k2,2n $bed_prox) <(sort -k1,1 -k2,2n $bed_dist) |
    groupBy -g 1,2,3,4,5,6,7,8 -c 13 -o sum |
    groupBy -g 1,2,3,4,5,6,7 -c 9 -o max > $tmp

intersectBed -sorted -wao -a $bed_te -b <(sort -k1,1 -k2,2n $bed_prox) <(sort -k1,1 -k2,2n $bed_dist) |
    groupBy -g 1,2,3,4,5,6,7,8 -c 13 -o sum |
    intersectBed -sorted -wao -f 1 -F 1 -a stdin -b $tmp | awk '$9==$(NF-1)' |
    groupBy -g 1,2,3,4,5,6,7 -c 8 -o first |
    awk 'BEGIN{OFS="\t"} $8==1 {$8="Proximal"} $8==2 {$8="Distal"} {print $1,$2,$3,$4,$5,$6,$7,$8}' > $out

