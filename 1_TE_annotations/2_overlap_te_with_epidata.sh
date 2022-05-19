#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=9-14
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
workdir=/scratch/twlab/hlee/zf_te
ID=$SLURM_ARRAY_TASK_ID

# TE classes
dir_feature=/scratch/twlab/hlee/genomes/danRer10
list=${dir_feature}/features.txt
feature=$( cat $list | sed "${ID}q;d" )
bed_te=${dir_feature}/${feature}.bed.gz
te=${feature#rmsk/danRer10.}


# INPUT DATA
tissue_list=${workdir}/tissues-e.txt
dir_in=${workdir}/0_cres


# OUTPUT
dir_out=1_annot
out=${dir_out}/epigenome_states_overlap_${te}.txt


# COMMANDS
echo -e "$te\toverlap\ttotal" > $out

for tissue_id in {1..11}
do
        tissue=$( cat $tissue_list | sed "${tissue_id}q;d" )
        out_tissue=${dir_out}/epigenome_states_overlap_${te}_in_${tissue}.txt

        prom=${dir_in}/activepromoter/${tissue}.active.promoter
        weak=${dir_in}/weakpromoter/${tissue}.weak.promoter
        enhc=${dir_in}/enhancer/${tissue}.active.enhancer
        hetero=${dir_in}/heterochromatin/${tissue}.Hetero-chromatin.bed

        prox=${dir_in}/proximalatacseq/${tissue}.proximal.open.bed
        dist=${dir_in}/distalatacseq/${tissue}.distal.open.bed

        umr=${dir_in}/umr/fylab_WGBS_zt_${tissue}_UMR.bed.gz
        lmr=${dir_in}/lmr/fylab_WGBS_zt_${tissue}_LMR.bed.gz

## Promoter / Weak promoter / Enhancer / Heterochromatin

	ovlp_prom=$( zcat $bed_te | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    intersectBed -a $prom -b stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	totl_prom=$( zcat $prom | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	echo -e "$ovlp_prom\t$totl_prom" > $out_tissue

	ovlp_enhc=$( zcat $bed_te | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    intersectBed -a $enhc -b stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	totl_enhc=$( zcat $enhc | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	echo -e "$ovlp_enhc\t$totl_enhc" >> $out_tissue

	ovlp_weak=$( zcat $bed_te | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    intersectBed -a $weak -b stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	totl_weak=$( zcat $weak | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	echo -e "$ovlp_weak\t$totl_weak" >> $out_tissue

	ovlp_hetero=$( zcat $bed_te | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    intersectBed -a $hetero -b stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	totl_hetero=$( zcat $hetero | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	echo -e "$ovlp_hetero\t$totl_hetero" >> $out_tissue

## ATAC peaks

	ovlp_prox=$( zcat $bed_te | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    intersectBed -a $prox -b stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	totl_prox=$( zcat $prox | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	echo -e "$ovlp_prox\t$totl_prox" >> $out_tissue

	ovlp_dist=$( zcat $bed_te | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    intersectBed -a $dist -b stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	totl_dist=$( zcat $dist | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	echo -e "$ovlp_dist\t$totl_dist" >> $out_tissue

## UMR / LMR

	ovlp_umr=$( zcat $bed_te | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    intersectBed -a $umr -b stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	totl_umr=$( zcat $umr | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	echo -e "$ovlp_umr\t$totl_umr" >> $out_tissue

	ovlp_lmr=$( zcat $bed_te | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    intersectBed -a $lmr -b stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	totl_lmr=$( zcat $lmr | sort -k1,1 -k2,2n | mergeBed -i stdin |
	    awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	echo -e "$ovlp_lmr\t$totl_lmr" >> $out_tissue

done

#out_list=$( cat $tissue_list | tr "\n" "," | sed "s/,$//g" )

paste ${dir_out}/epigenome_states_overlap_${te}_in_{blood,brain,colon,heart,intestine,kidney,liver,muscle,skin,spleen,testis}.txt |
        awk -F"\t" '{te=0; sum=0; for (i=1;i<NF;i+=2) {te+=$i; sum+=$(i+1); percnt=sprintf("%.3f", te/sum)} print percnt"\t"te"\t"sum}' >> $out


