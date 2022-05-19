#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=1-11
#SBATCH --job-name=propor_ind_te

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID


# TE 6 bed files
te=danRer10.TE.bed.gz
dna=danRer10.DNA.bed.gz
ltr=danRer10.LTR.bed.gz
line=danRer10.LINE.bed.gz
sine=danRer10.SINE.bed.gz
rc=danRer10.RC.bed.gz

te_bed=$(ls danRer10.{TE,DNA,LTR,LINE,SINE,RC}_frag.bed.gz )


# INPUT DATA
list=tissues-e.txt
tissue=$( cat $list | sed "${ID}q;d" )

dir_in=0_cres
prom=${dir_in}/activepromoter/${tissue}.active.promoter
weak=${dir_in}/weakpromoter/${tissue}.weak.promoter
enhc=${dir_in}/enhancer/${tissue}.active.enhancer
hetero=${dir_in}/heterochromatin/${tissue}.Hetero-chromatin.bed

prox=${dir_in}/proximalatacseq/${tissue}.proximal.open.bed
dist=${dir_in}/distalatacseq/${tissue}.distal.open.bed

umr=${dir_in}/umr/fylab_WGBS_zt_${tissue}_UMR.bed.gz
lmr=${dir_in}/lmr/fylab_WGBS_zt_${tissue}_LMR.bed.gz


# OUTPUT
dir_out=2_dynamic
out=${dir_out}/fraction_epigenome_states_of_TE_in_${tissue}.txt


# COMMANDS

for bed in $te_bed
do

	name=${bed##*/}
	name=${name%.bed.gz}
	name=${name#danRer10.}

	cnt_te_prom=$( intersectBed -u -a $bed -b $prom | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_prom=$( awk -v overlap="$cnt_te_prom" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_prom}\t${cnt_te_prom}\t${cnt_te}\tActive_promoter\t${name}\t$tissue"

	cnt_te_weak=$( intersectBed -u -a $bed -b $weak | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_weak=$( awk -v overlap="$cnt_te_weak" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_weak}\t${cnt_te_weak}\t${cnt_te}\tWeak_promoter\t${name}\t$tissue"

	cnt_te_enhc=$( intersectBed -u -a $bed -b $enhc | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_enhc=$( awk -v overlap="$cnt_te_enhc" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_enhc}\t${cnt_te_enhc}\t${cnt_te}\tActive_enhancer\t${name}\t$tissue"

	cnt_te_hetero=$( intersectBed -u -a $bed -b $hetero | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_hetero=$( awk -v overlap="$cnt_te_hetero" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_hetero}\t${cnt_te_hetero}\t${cnt_te}\tHeterochromatin\t${name}\t$tissue"

	cnt_te_prox=$( intersectBed -u -a $bed -b $prox | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_prox=$( awk -v overlap="$cnt_te_prox" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_prox}\t${cnt_te_prox}\t${cnt_te}\tProximal_ATAC\t${name}\t$tissue"

	cnt_te_dist=$( intersectBed -u -a $bed -b $dist | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_dist=$( awk -v overlap="$cnt_te_dist" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_dist}\t${cnt_te_dist}\t${cnt_te}\tDistal_ATAC\t${name}\t$tissue"

	cnt_te_umr=$( intersectBed -u -a $bed -b $umr | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_umr=$( awk -v overlap="$cnt_te_umr" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_umr}\t${cnt_te_umr}\t${cnt_te}\tUMR\t${name}\t$tissue"

	cnt_te_lmr=$( intersectBed -u -a $bed -b $lmr | wc -l )
	cnt_te=$( zcat $bed | wc -l )
	rate_lmr=$( awk -v overlap="$cnt_te_lmr" -v total="$cnt_te" 'BEGIN{ printf("%.5f", overlap/total) }' )
	echo -e "${rate_lmr}\t${cnt_te_lmr}\t${cnt_te}\tLMR\t${name}\t$tissue"

done > $out
