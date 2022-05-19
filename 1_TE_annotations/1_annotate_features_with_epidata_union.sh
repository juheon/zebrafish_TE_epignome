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


# OUTPUT
dir_out=1_annot
feature=${feature##*/}
feature=${feature#danRer10.}
feature=${feature#GRCz10.91.}
out=${dir_out}/${feature}.annotation_epigenome_union.txt


# COMMANDS
echo $feature > $out


prom=${dir_in}/active_promoters_union.bed.gz
weak=${dir_in}/weak_promoters_union.bed.gz
enhc=${dir_in}/active_enhancers_union.bed.gz
hetero=${dir_in}/heterochromatin_union.bed.gz

prox=${dir_in}/proximal_atac_union.bed.gz
dist=${dir_in}/distal_atac_union.bed.gz

umr=${dir_in}/UMRs_union.bed.gz
lmr=${dir_in}/LMRs_union.bed.gz

## Promoter, Enhancer, Weak promoter, Heterochromatin
intersectBed -a $bed_feature -b $prom | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out

intersectBed -a $bed_feature -b $enhc | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out

intersectBed -a $bed_feature -b $weak | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out

intersectBed -a $bed_feature -b $hetero | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out

## ATAC peaks
intersectBed -a $bed_feature -b $prox | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out

intersectBed -a $bed_feature -b $dist | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out

## UMR / LMR
intersectBed -a $bed_feature -b $umr | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out

intersectBed -a $bed_feature -b $lmr | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out

## Genomic feature size
zcat $bed_feature | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out


