#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=1-11
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID

# GENOME FEATURES
dir_te=/scratch/twlab/hlee/genomes/danRer10/rmsk
te=${dir_te}/danRer10.TE.bed.gz
dna=${dir_te}/danRer10.DNA.bed.gz
ltr=${dir_te}/danRer10.LTR.bed.gz
line=${dir_te}/danRer10.LINE.bed.gz
sine=${dir_te}/danRer10.SINE.bed.gz
rc=${dir_te}/danRer10.RC.bed.gz

te_bed=$(ls ${dir_te}/danRer10.{TE,DNA,LTR,LINE,SINE,RC}.bed.gz )


# INPUT DATA
list=tissues-e.txt
tissue=$( cat $list | sed "${ID}q;d" )

dir_in=/scratch/twlab/hlee/zf_te/0_cres
prom=${dir_in}/activepromoter/${tissue}.active.promoter
weak=${dir_in}/weakpromoter/${tissue}.weak.promoter
enhc=${dir_in}/enhancer/${tissue}.active.enhancer
hetero=${dir_in}/heterochromatin/${tissue}.Hetero-chromatin.bed

prox=${dir_in}/proximalatacseq/${tissue}.proximal.open.bed
dist=${dir_in}/distalatacseq/${tissue}.distal.open.bed

umr=${dir_in}/umr/fylab_WGBS_zt_${tissue}_UMR.bed.gz
lmr=${dir_in}/lmr/fylab_WGBS_zt_${tissue}_LMR.bed.gz


# OUTPUT
dir_out=1_annot
out=${dir_out}/TE_contribution_epi_states_${tissue}.txt



# COMMANDS
#echo "$tissue" > $out

for bed in $te_bed
do

    name=${bed##*/}
    name=${name%.bed.gz}
    name=${name#danRer10.}

    cnt_prom_te=$( intersectBed -a $prom -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    cnt_prom=$( zcat $prom | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    rate_prom=$( awk -v overlap="$cnt_prom_te" -v total="$cnt_prom" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_prom}\t${cnt_prom_te}\t${cnt_prom}\tActive_promoter\t${name}\t$tissue"

    cnt_weak_te=$( intersectBed -sorted -a $weak -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    cnt_weak=$( zcat $weak | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    rate_weak=$( awk -v overlap="$cnt_weak_te" -v total="$cnt_weak" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_weak}\t${cnt_weak_te}\t${cnt_weak}\tWeak_promoter\t${name}\t$tissue"

    cnt_enhc_te=$( intersectBed -a $enhc -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    cnt_enhc=$( zcat $enhc | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    rate_enhc=$( awk -v overlap="$cnt_enhc_te" -v total="$cnt_enhc" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_enhc}\t${cnt_enhc_te}\t${cnt_enhc}\tActive_enhancer\t${name}\t$tissue"

    cnt_hetero_te=$( intersectBed -sorted -a $hetero -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    cnt_hetero=$( zcat $hetero | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    rate_hetero=$( awk -v overlap="$cnt_hetero_te" -v total="$cnt_hetero" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_hetero}\t${cnt_hetero_te}\t${cnt_hetero}\tHeterochromatin\t${name}\t$tissue"

    cnt_prox_te=$( intersectBed -a $prox -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    cnt_prox=$( zcat $prox | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    rate_prox=$( awk -v overlap="$cnt_prox_te" -v total="$cnt_prox" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_prox}\t${cnt_prox_te}\t${cnt_prox}\tProximal_ATAC\t${name}\t$tissue"

    cnt_dist_te=$( intersectBed -sorted -a $dist -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    cnt_dist=$( zcat $dist | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    rate_dist=$( awk -v overlap="$cnt_dist_te" -v total="$cnt_dist" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_dist}\t${cnt_dist_te}\t${cnt_dist}\tDistant_ATAC\t${name}\t$tissue"

    cnt_umr_te=$( intersectBed -sorted -a $umr -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    cnt_umr=$( zcat $umr | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    rate_umr=$( awk -v overlap="$cnt_umr_te" -v total="$cnt_umr" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_umr}\t${cnt_umr_te}\t${cnt_umr}\tUMR\t${name}\t$tissue"

    cnt_lmr_te=$( intersectBed -sorted -a $lmr -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    cnt_lmr=$( zcat $lmr | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
    rate_lmr=$( awk -v overlap="$cnt_lmr_te" -v total="$cnt_lmr" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_lmr}\t${cnt_lmr_te}\t${cnt_lmr}\tLMR\t${name}\t$tissue"

done > $out


