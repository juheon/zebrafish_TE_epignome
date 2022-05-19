#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=9-14
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID

# TE
dir_te=/scratch/twlab/hlee/genomes/danRer10
list=${dir_te}/features.txt

te=$( cat $list | sed "${ID}q;d" )
te=${te/rmsk\/danRer10./}

bed=${dir_te}/rmsk/danRer10.${te}.bed.gz


# INPUT DATA
dir_in=/scratch/twlab/hlee/zf_te/0_cres

prom=${dir_in}/active_promoters_union.bed.gz
weak=${dir_in}/weak_promoters_union.bed.gz
enhc=${dir_in}/active_enhancers_union.bed.gz
hetero=${dir_in}/heterochromatin_union.bed.gz

prox=${dir_in}/proximal_atac_union.bed.gz
dist=${dir_in}/distal_atac_union.bed.gz

umr=${dir_in}/UMRs_union.bed.gz
lmr=${dir_in}/LMRs_union.bed.gz


# OUTPUT
dir_out=1_annot
out=${dir_out}/TE_contribution_${te}_epi_states_union.txt


# COMMANDS
echo "$te" > $out

cnt_prom_te=$( intersectBed -a $prom -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
cnt_prom=$( zcat $prom | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
rate_prom=$( awk -v overlap="$cnt_prom_te" -v total="$cnt_prom" 'BEGIN{ printf("%.4f", overlap/total) }' )
echo -e "${rate_prom}\t${cnt_prom_te}\t${cnt_prom}\tActive_promoter" >> $out

cnt_weak_te=$( intersectBed -sorted -a $weak -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
cnt_weak=$( zcat $weak | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
rate_weak=$( awk -v overlap="$cnt_weak_te" -v total="$cnt_weak" 'BEGIN{ printf("%.4f", overlap/total) }' )
echo -e "${rate_weak}\t${cnt_weak_te}\t${cnt_weak}\tWeak_promoter" >> $out

cnt_enhc_te=$( intersectBed -a $enhc -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
cnt_enhc=$( zcat $enhc | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
rate_enhc=$( awk -v overlap="$cnt_enhc_te" -v total="$cnt_enhc" 'BEGIN{ printf("%.4f", overlap/total) }' )
echo -e "${rate_enhc}\t${cnt_enhc_te}\t${cnt_enhc}\tActive_enhancer" >> $out

cnt_hetero_te=$( intersectBed -sorted -a $hetero -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
cnt_hetero=$( zcat $hetero | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
rate_hetero=$( awk -v overlap="$cnt_hetero_te" -v total="$cnt_hetero" 'BEGIN{ printf("%.4f", overlap/total) }' )
echo -e "${rate_hetero}\t${cnt_hetero_te}\t${cnt_hetero}\tHeterochromatin" >> $out

cnt_prox_te=$( intersectBed -a $prox -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
cnt_prox=$( zcat $prox | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
rate_prox=$( awk -v overlap="$cnt_prox_te" -v total="$cnt_prox" 'BEGIN{ printf("%.4f", overlap/total) }' )
echo -e "${rate_prox}\t${cnt_prox_te}\t${cnt_prox}\tProximal_ATAC" >> $out

cnt_dist_te=$( intersectBed -sorted -a $dist -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
cnt_dist=$( zcat $dist | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
rate_dist=$( awk -v overlap="$cnt_dist_te" -v total="$cnt_dist" 'BEGIN{ printf("%.4f", overlap/total) }' )
echo -e "${rate_dist}\t${cnt_dist_te}\t${cnt_dist}\tDistant_ATAC" >> $out

cnt_umr_te=$( intersectBed -sorted -a $umr -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
cnt_umr=$( zcat $umr | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
rate_umr=$( awk -v overlap="$cnt_umr_te" -v total="$cnt_umr" 'BEGIN{ printf("%.4f", overlap/total) }' )
echo -e "${rate_umr}\t${cnt_umr_te}\t${cnt_umr}\tUMR" >> $out

cnt_lmr_te=$( intersectBed -sorted -a $lmr -b $bed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
cnt_lmr=$( zcat $lmr | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
rate_lmr=$( awk -v overlap="$cnt_lmr_te" -v total="$cnt_lmr" 'BEGIN{ printf("%.4f", overlap/total) }' )
echo -e "${rate_lmr}\t${cnt_lmr_te}\t${cnt_lmr}\tLMR" >> $out


