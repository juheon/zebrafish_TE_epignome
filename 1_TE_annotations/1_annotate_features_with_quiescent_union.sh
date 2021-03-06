#!/bin/bash
# Author: Hyung Joo Lee

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID
workdir=`pwd`

# GENOME FEATURES
list=features.txt
feature=$( cat $list | sed "${ID}q;d" )
bed_feature=${feature}.bed.gz

# INPUT DATA
dir_in=${workdir}/0_cres/union

# OUTPUT
dir_out=1_annot
feature=${feature##*/}
feature=${feature#danRer10.}
feature=${feature#GRCz10.91.}
out=${dir_out}/${feature}.annotation_quiescent_union.txt

# COMMANDS
echo $feature > $out

quie=${dir_in}/quiescent_union.bed.gz

intersectBed -a $bed_feature -b $quie | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out

## Genomic feature size
zcat $bed_feature | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out


