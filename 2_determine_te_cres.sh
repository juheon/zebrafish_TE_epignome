#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --array=1-11
#SBATCH --mem=2G
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=te_cres_tissue

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
bed_prom=${dir_in}/activepromoter/${tissue}.active.promoter
bed_weak=${dir_in}/weakpromoter/${tissue}.weak.promoter
bed_enhc=${dir_in}/enhancer/${tissue}.active.enhancer
bed_hetero=${dir_in}/heterochromatin/${tissue}.Hetero-chromatin.bed


# OUTPUT
dir_out=2_dynamic
out=${dir_out}/danRer10.TEs_${tissue}_CREs.bed.gz


# COMMANDS
intersectBed -sorted -wao -a $bed_te -b <(sort -k1,1 -k2,2n $bed_prom) <(sort -k1,1 -k2,2n $bed_weak) <(sort -k1,1 -k2,2n $bed_enhc) <(sort -k1,1 -k2,2n $bed_hetero) |
    groupBy -g 1,2,3,8 -c 12 -o sum |
    groupBy -c 5 -o max -full |
    awk 'BEGIN{OFS="\t"} $4==1 {$4="Active_promoter"} $4==2 {$4="Weak_promoter"}
			 $4==3 {$4="Active_enhancer"} $4==4 {$4="Heterochromatin"} {print $1,$2,$3,$4}' | gzip -nc > $out

