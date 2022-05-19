#!/bin/bash
# Author: Hyung Joo Lee

##SBATCH --mem=2G
#SBATCH --array=1-14
#SBATCH --job-name=annot_feat_hic

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID
workdir=`pwd`

# GENOME FEATURES
features.txt
feature=$( cat $list | sed "${ID}q;d" )
bed_feature=${feature}.bed.gz


# INPUT DATA
tissue_list=${workdir}/tissues-e.txt
dir_in=${workdir}/0_cres


# OUTPUT
dir_out=${workdir}/1_annot
feature=${feature##*/}
feature=${feature#danRer10.}
feature=${feature#GRCz10.91.}
#out=${dir_out}/${feature}.annotation_ctcf_atac.txt
out=${dir_out}/${feature}.annotation_atac_ctcf.txt


# COMMANDS
echo $feature > $out

for tissue_id in {2,8}
do
	tissue=$( cat $tissue_list | sed "${tissue_id}q;d" )
#	out_tissue=${dir_out}/${feature}_${tissue}.annotation_ctcf_atac.txt

#	tad=${dir_in}/TADs/TAD_boundaries_${tissue}.ctcf_atac.bed.gz
#	anchor=${dir_in}/loop_anchors/loop_anchors_${tissue}.ctcf_atac.bed.gz

	out_tissue=${dir_out}/${feature}_${tissue}.annotation_atac_ctcf.txt

	tad=${dir_in}/TADs/TAD_boundaries_${tissue}.atac_ctcf.bed.gz
	anchor=${dir_in}/loop_anchors/loop_anchors_${tissue}.atac_ctcf.bed.gz

## TAD boundaris and loop anchors 
	intersectBed -a $bed_feature -b $tad | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' > $out_tissue

	intersectBed -a $bed_feature -b $anchor | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out_tissue

## Genomic feature size
	zcat $bed_feature | sort -k1,1 -k2,2n |
	    mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' >> $out_tissue

done


#out_list=$( cat $tissue_list | tr "\n" "," | sed "s/^/{/g; s/,$/}/g" )
#out_list=${dir_out}/${feature}_${out_list}.annotation_epigenome.txt

#paste ${dir_out}/${feature}_{brain,muscle}.annotation_ctcf_atac.txt |
paste ${dir_out}/${feature}_{brain,muscle}.annotation_atac_ctcf.txt |
	awk -F"\t" '{s=0; for (i=1;i<=NF;i++) {s+=$i} print s}' >> $out

