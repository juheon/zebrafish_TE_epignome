#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --array=2,8
#SBATCH --mem=2G
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=te_cres_tissue

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID

workdir=/scratch/twlab/hlee/zf_te

# TE classes
dir_te=/scratch/twlab/hlee/genomes/danRer10/rmsk
bed_te=${dir_te}/danRer10.TE_frag.bed.gz


# INPUT DATA
list=${workdir}/tissues-e.txt
tissue=$( cat $list | sed "${ID}q;d" )

dir_in=${workdir}/0_cres
bed_boundary=${dir_in}/TADs/TAD_boundaries_${tissue}.ctcf_atac.bed.gz
bed_anchor=${dir_in}/loop_anchors/loop_anchors_${tissue}.ctcf_atac.bed.gz


# OUTPUT
dir_out=${workdir}/2_dynamic
out1=${dir_out}/danRer10.TEs_${tissue}_boundaries.ctcf_atac.bed.gz
out2=${dir_out}/danRer10.TEs_${tissue}_anchors.ctcf_atac.bed.gz


# COMMANDS
intersectBed -sorted -wao -a $bed_te -b $bed_boundary |
    awk 'BEGIN{OFS="\t"} $14>0 {$8="TAD_boundary"} $14==0 {$8="No"} {print $1,$2,$3,$4,$5,$6,$7,$8}' |
	groupBy -g 1,2,3,4,5,6,7 -c 8 -o distinct |
	awk -F"\t" -vOFS="\t" '$8~"," {$8="TAD_boundary"} {print $0}' | gzip -nc > $out1

intersectBed -sorted -wao -a $bed_te -b $bed_anchor |
    awk 'BEGIN{OFS="\t"} $14>0 {$8="loop_anchor"} $14==0 {$8="No"} {print $1,$2,$3,$4,$5,$6,$7,$8}' |
	groupBy -g 1,2,3,4,5,6,7 -c 8 -o distinct |  
        awk -F"\t" -vOFS="\t" '$8~"," {$8="loop_anchor"} {print $0}' | gzip -nc > $out2

