#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID


# TE classes
dir_te=/scratch/twlab/hlee/genomes/danRer10/rmsk
bed_te=${dir_te}/danRer10.TEs.bed.gz

# INPUT DATA
cre=2_dynamic/danRer10.TEs_11_tissues_CREs.bed.gz
atac=2_dynamic/danRer10.TEs_11_tissues_ATAC.bed.gz
umr_lmr=2_dynamic/danRer10.TEs_11_tissues_UMR_LMR.bed.gz


# OUTPUT
dir_out=2_dynamic
out_cre=${dir_out}/cnt_te_cre_shared_tissue.txt
out_atac=${dir_out}/cnt_te_atac_shared_tissue.txt
out_umr_lmr=${dir_out}/cnt_te_umr_lmr_shared_tissue.txt


# COMMANDS
zcat $cre |
    awk -F"\t" 'BEGIN{ for (i=0; i<=11; i++) {cnt_prom[i]=0; cnt_weak[i]=0; cnt_enhc[i]=0; cnt_hete[i]=0 } }
        { prom=0; weak=0; enhc=0; hete=0;
        for(i=4; i<=NF; i++) {
            if( $i=="Active_promoter" ) {prom++}
            if( $i=="Weak_promoter" ) {weak++}
            if( $i=="Active_enhancer" ) {enhc++}
            if( $i=="Heterochromatin" ) {hete++} }
        cnt_prom[prom]++; cnt_weak[weak]++; cnt_enhc[enhc]++; cnt_hete[hete]++ }
        END{ OFS="\t"; print "Active_promoter", "Weak_promoter", "Active_enhancer", "Heterochromatin";
        for(i=1; i<=11; i++) { print cnt_prom[i],cnt_weak[i],cnt_enhc[i],cnt_hete[i] } } ' > $out_cre

zcat $atac |
    awk -F"\t" 'BEGIN{ for (i=0; i<=11; i++) {cnt_prox[i]=0; cnt_dist[i]=0 } }
        { prox=0; dist=0;
        for(i=4; i<=NF; i++) {
            if( $i=="Proximal" ) {prox++}
            if( $i=="Distal" ) {dist++} }
        cnt_prox[prox]++; cnt_dist[dist]++ }
        END{ OFS="\t"; print "Proximal", "Distal";
        for(i=1; i<=11; i++) { print cnt_prox[i],cnt_dist[i] } } ' > $out_atac

zcat $umr_lmr |
    awk -F"\t" 'BEGIN{ for (i=0; i<=11; i++) {cnt_umr[i]=0; cnt_lmr[i]=0 } }
        { umr=0; lmr=0;
        for(i=4; i<=NF; i++) {
            if( $i=="UMR" ) {umr++}
            if( $i=="LMR" ) {lmr++} }
        cnt_umr[umr]++; cnt_lmr[lmr]++ }
        END{ OFS="\t"; print "UMR", "LMR";
        for(i=1; i<=11; i++) { print cnt_umr[i],cnt_lmr[i] } } ' > $out_umr_lmr



