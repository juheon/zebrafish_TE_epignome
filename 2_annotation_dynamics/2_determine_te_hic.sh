#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --array=2,8
#SBATCH --mem=2G
#SBATCH --job-name=te_cres_tissue

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID

workdir=`pwd`

# TE classes
bed_te=danRer10.TE_frag.bed.gz


# INPUT DATA
list=${workdir}/tissues-e.txt
tissue=$( cat $list | sed "${ID}q;d" )

dir_in=${workdir}/0_cres
bed_boundary=${dir_in}/TADs/TAD_boundaries_${tissue}.bed.gz
bed_anchor=${dir_in}/loop_anchors/loop_anchors_${tissue}.bed.gz


# OUTPUT
dir_out=${workdir}/2_dynamic
out1=${dir_out}/danRer10.TEs_${tissue}_boundaries.bed.gz
out2=${dir_out}/danRer10.TEs_${tissue}_anchors.bed.gz


# COMMANDS
intersectBed -sorted -wao -a $bed_te -b $bed_boundary |
    awk 'BEGIN{OFS="\t"} $11>0 {$4="TAD_boundary"} $11==0 {$4="No"} {print $1,$2,$3,$4}' | gzip -nc > $out1

intersectBed -sorted -wao -a $bed_te -b $bed_anchor |
    awk 'BEGIN{OFS="\t"} $11>0 {$4="loop_anchor"} $11==0 {$4="No"} {print $1,$2,$3,$4}' | gzip -nc > $out2

