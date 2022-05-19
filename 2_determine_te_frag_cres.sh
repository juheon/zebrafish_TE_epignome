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
tmp=${dir_out}/danRer10.TEs_${tissue}_CREs.tmp.bed
out=${dir_out}/danRer10.TEs_${tissue}_CREs.bed


# COMMANDS
intersectBed -sorted -wao -a $bed_te -b <(sort -k1,1 -k2,2n $bed_prom) <(sort -k1,1 -k2,2n $bed_weak) <(sort -k1,1 -k2,2n $bed_enhc) <(sort -k1,1 -k2,2n $bed_hetero) |
    groupBy -g 1,2,3,4,5,6,7,8 -c 12 -o sum |
    groupBy -g 1,2,3,4,5,6,7 -c 9 -o max > $tmp

intersectBed -sorted -wao -a $bed_te -b <(sort -k1,1 -k2,2n $bed_prom) <(sort -k1,1 -k2,2n $bed_weak) <(sort -k1,1 -k2,2n $bed_enhc) <(sort -k1,1 -k2,2n $bed_hetero) |
    groupBy -g 1,2,3,4,5,6,7,8 -c 12 -o sum |
    intersectBed -sorted -wao -f 1 -F 1 -a stdin -b $tmp | awk '$9==$(NF-1)' |
    groupBy -g 1,2,3,4,5,6,7 -c 8 -o first |
    awk 'BEGIN{OFS="\t"} $8==1 {$8="Active_promoter"} $8==2 {$8="Weak_promoter"}
			 $8==3 {$8="Active_enhancer"} $8==4 {$8="Heterochromatin"} {print $1,$2,$3,$4,$5,$6,$7,$8}' | gzip -nc > $out

#rm $tmp
