#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=9-14
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
workdir=`pwd`
ID=$SLURM_ARRAY_TASK_ID

# TE classes
list=TEfeatures.txt
feature=$( cat $list | sed "${ID}q;d" )
bed_te=${feature}.bed.gz
te=${feature#rmsk/danRer10.}


# INPUT DATA
tissue_list=${workdir}/tissues-e.txt
dir_in=${workdir}/0_cres


# OUTPUT
dir_out=${workdir}/1_annot
#out=${dir_out}/queiscent_states_overlap_${te}.txt
out=${dir_out}/atac_ctcf_zf_states_overlap_${te}.txt


# COMMANDS
echo -e "$te\toverlap\ttotal" > $out

for tissue_id in {2,8}
do
        tissue=$( cat $tissue_list | sed "${tissue_id}q;d" )
        out_tissue=${dir_out}/atac_ctcf_zf_states_overlap_${te}_in_${tissue}.txt
	bound=${dir_in}/TADs/TAD_boundaries_${tissue}.atac_ctcf_zebrafish.bed.gz
	anchors=${dir_in}/loop_anchors/loop_anchors_${tissue}.cres_ctcf_zebrafish.bed.gz

## TAD boundaries / Loop anchors

	ovlp_bound=$( zcat $bed_te | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    intersectBed -a $bound -b stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	totl_bound=$( zcat $bound | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	echo -e "$ovlp_bound\t$totl_bound" >> $out_tissue

	ovlp_anchors=$( zcat $bed_te | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    intersectBed -a $anchors -b stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	totl_anchors=$( zcat $anchors | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	echo -e "$ovlp_anchors\t$totl_anchors" >> $out_tissue



done

paste ${dir_out}/atac_ctcf_zf_states_overlap_${te}_in_{brain,muscle}.txt |
        awk -F"\t" '{te=0; sum=0; for (i=1;i<NF;i+=2) {te+=$i; sum+=$(i+1); percnt=sprintf("%.3f", te/sum)} print percnt"\t"te"\t"sum}' >> $out



