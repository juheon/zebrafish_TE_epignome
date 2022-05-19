#!/bin/bash
# Author: Hyung Joo Lee

##SBATCH --mem=2G
#SBATCH --array=14
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
tissue_list=${workdir}/tissues-e.txt
dir_in=${workdir}/0_cres


# OUTPUT
dir_out=1_annot
feature=${feature##*/}
feature=${feature#danRer10.}
feature=${feature#GRCz10.91.}
out=${dir_out}/${feature}.annotation_epigenome.txt


# COMMANDS
echo $feature > $out

for tissue_id in {1..11}
do
	tissue=$( cat $tissue_list | sed "${tissue_id}q;d" )
	out_tissue=${dir_out}/${feature}_${tissue}.annotation_epigenome.txt

	prom=${dir_in}/activepromoter/${tissue}.active.promoter
	weak=${dir_in}/weakpromoter/${tissue}.weak.promoter
	enhc=${dir_in}/enhancer/${tissue}.active.enhancer
	hetero=${dir_in}/heterochromatin/${tissue}.Hetero-chromatin.bed

	prox=${dir_in}/proximalatacseq/${tissue}.proximal.open.bed
	dist=${dir_in}/distalatacseq/${tissue}.distal.open.bed

	umr=${dir_in}/umr/fylab_WGBS_zt_${tissue}_UMR.bed.gz
	lmr=${dir_in}/lmr/fylab_WGBS_zt_${tissue}_LMR.bed.gz

## Promoter -> Enhancer -> Weak promoter -> Enhancer -> Heterochromatin
	intersectBed -a $bed_feature -b $prom | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' > $out_tissue

	intersectBed -a $bed_feature -b $enhc |
	    intersectBed -v -a stdin -b $prom | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out_tissue

	intersectBed -a $bed_feature -b $weak |
            intersectBed -v -a stdin -b $prom $enhc | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out_tissue

	intersectBed -a $bed_feature -b $hetero |
            intersectBed -v -a stdin -b $prom $enhc $weak | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out_tissue

## ATAC peaks
	intersectBed -a $bed_feature -b $prox | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out_tissue


	intersectBed -a $bed_feature -b $dist | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out_tissue

## UMR / LMR
	intersectBed -a $bed_feature -b $umr | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out_tissue

	intersectBed -a $bed_feature -b $lmr | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out_tissue

## Genomic feature size
	zcat $bed_feature | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out_tissue

done


#out_list=$( cat $tissue_list | tr "\n" "," | sed "s/^/{/g; s/,$/}/g" )
#out_list=${dir_out}/${feature}_${out_list}.annotation_epigenome.txt

paste ${dir_out}/${feature}_{blood,brain,colon,heart,intestine,kidney,liver,muscle,skin,spleen,testis}.annotation_epigenome.txt |
	awk -F"\t" '{s=0; for (i=1;i<=NF;i++) {s+=$i} print s}' >> $out

