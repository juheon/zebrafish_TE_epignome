#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=9-14
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID

# TE
dir_te=/scratch/twlab/hlee/genomes/danRer10
list=${dir_te}/features.txt

te=$( cat $list | sed "${ID}q;d" )
te=${te/rmsk\/danRer10./}

bed=${dir_te}/rmsk/danRer10.${te}.bed.gz


# INPUT DATA
dir_in=/scratch/twlab/hlee/zf_te/0_cres

quie=${dir_in}/quiescent_union.bed.gz


# OUTPUT
dir_out=1_annot
out=${dir_out}/TE_contribution_${te}_quiescent_union.txt


# COMMANDS
echo "$te" > $out

cnt_quie_te=$( intersectBed -sorted -a $quie -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
cnt_quie=$( zcat $quie | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
rate_quie=$( awk -v overlap="$cnt_quie_te" -v total="$cnt_quie" 'BEGIN{ printf("%.4f", overlap/total) }' )
echo -e "${rate_quie}\t${cnt_quie_te}\t${cnt_quie}\tQuiescent" >> $out


