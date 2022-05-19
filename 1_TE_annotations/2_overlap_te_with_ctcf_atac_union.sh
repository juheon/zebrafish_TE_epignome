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
dir_in=/scratch/twlab/hlee/zf_te/0_cres/union
boundaries=${dir_in}/TAD_boundaries_union.ctcf_atac.bed.gz
anchors=${dir_in}/loop_anchors_union.ctcf_atac.bed.gz


# OUTPUT
dir_out=1_annot
out=${dir_out}/TE_contribution_${te}_ctcf_atac_union.txt


# COMMANDS
echo "$te" > $out

cnt_boundaries_te=$( intersectBed -sorted -a $boundaries -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
cnt_boundaries=$( zcat $boundaries | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
rate_boundaries=$( awk -v overlap="$cnt_boundaries_te" -v total="$cnt_boundaries" 'BEGIN{ printf("%.4f", overlap/total) }' )
echo -e "${rate_boundaries}\t${cnt_boundaries_te}\t${cnt_boundaries}\tUMR" >> $out

cnt_anchors_te=$( intersectBed -sorted -a $anchors -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
cnt_anchors=$( zcat $anchors | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
rate_anchors=$( awk -v overlap="$cnt_anchors_te" -v total="$cnt_anchors" 'BEGIN{ printf("%.4f", overlap/total) }' )
echo -e "${rate_anchors}\t${cnt_anchors_te}\t${cnt_anchors}\tLMR" >> $out


