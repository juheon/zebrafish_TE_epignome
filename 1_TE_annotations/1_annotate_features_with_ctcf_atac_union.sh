#!/bin/bash
# Author: Hyung Joo Lee

##SBATCH --mem=2G
#SBATCH --array=1-14
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID
workdir=`pwd`

# GENOME FEATURES
list=features.txt
feature=$( cat $list | sed "${ID}q;d" )
bed_feature=${feature}.bed.gz

# INPUT DATA
dir_in=${workdir}/0_cres

#boundaries=${dir_in}/union/TAD_boundaries_union.ctcf_atac.bed.gz
#anchors=${dir_in}/union/loop_anchors_union.ctcf_atac.bed.gz
boundaries=${dir_in}/union/TAD_boundaries_union.atac_ctcf.bed.gz
anchors=${dir_in}/union/loop_anchors_union.atac_ctcf.bed.gz


# OUTPUT
dir_out=${workdir}/1_annot
feature=${feature##*/}
feature=${feature#danRer10.}
feature=${feature#GRCz10.91.}
#out=${dir_out}/${feature}.annotation_ctcf_atac_union.txt
out=${dir_out}/${feature}.annotation_atac_ctcf_union.txt


# COMMANDS
echo $feature > $out

## TAD boundaries and loop anchors
intersectBed -a $bed_feature -b $boundaries | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out

intersectBed -a $bed_feature -b $anchors | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out

## Genomic feature size
zcat $bed_feature | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out


