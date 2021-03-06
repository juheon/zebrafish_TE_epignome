#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=9-14
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID

# TE classes
list=TEfeatures.txt
feature=$( cat $list | sed "${ID}q;d" )
bed_te=${dir_feature}/${feature}.bed.gz
te=${feature#rmsk/danRer10.}


# INPUT DATA
tissue_list=tissues-e.txt
dir_in=0_cres


# OUTPUT
dir_out=1_annot
#out=${dir_out}/queiscent_states_overlap_${te}.txt
out=${dir_out}/hic_states_overlap_${te}.txt


# COMMANDS
echo -e "$te\toverlap\ttotal" > $out

for tissue_id in {2,8}
do
        tissue=$( cat $tissue_list | sed "${tissue_id}q;d" )
        out_tissue=${dir_out}/hic_states_overlap_${te}_in_${tissue}.txt

	quiescent=${dir_in}/quiescent/quiescent_${tissue}.bed.gz
	bound=${dir_in}/TADs/TAD_boundaries_${tissue}.bed.gz
	anchors=${dir_in}/loop_anchors/loop_anchors_${tissue}.bed.gz


## Promoter / Weak promoter / Enhancer / Heterochromatin

#	ovlp_quies=$( zcat $bed_te | sort -k1,1 -k2,2n | mergeBed -i stdin |
#	    intersectBed -a $quiescent -b stdin |
#	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
#	totl_quies=$( zcat $quiescent | sort -k1,1 -k2,2n | mergeBed -i stdin |
#	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
#	echo -e "$ovlp_quies\t$totl_quies" >> $out_tissue

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

#out_list=$( cat $tissue_list | tr "\n" "," | sed "s/,$//g" )

#paste ${dir_out}/quiescent_states_overlap_${te}_in_{blood,brain,colon,heart,intestine,kidney,liver,muscle,skin,spleen,testis}.txt |
#        awk -F"\t" '{te=0; sum=0; for (i=1;i<NF;i+=2) {te+=$i; sum+=$(i+1); percnt=sprintf("%.3f", te/sum)} print percnt"\t"te"\t"sum}' >> $out

paste ${dir_out}/hic_states_overlap_${te}_in_{brain,muscle}.txt |
        awk -F"\t" '{te=0; sum=0; for (i=1;i<NF;i+=2) {te+=$i; sum+=$(i+1); percnt=sprintf("%.3f", te/sum)} print percnt"\t"te"\t"sum}' >> $out



