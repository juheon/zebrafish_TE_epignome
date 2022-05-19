#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=te_cres_tissue

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID


# TE classes
dir_te=/scratch/twlab/hlee/genomes/danRer10/rmsk
bed_te=${dir_te}/danRer10.TE_frag.bed.gz


# INPUT DATA
atac=2_dynamic/danRer10.TEs_*_ATAC.bed
cre=2_dynamic/danRer10.TEs_*_CREs.bed
umr_lmr=2_dynamic/danRer10.TEs_*_UMR_LMR.bed
dname=2_dynamic/danRer10.TE_wCG.DNAme_*_cov5.bed.gz
#boundary=2_dynamic/danRer10.TEs_*_boundaries.bed.gz
#anchor=2_dynamic/danRer10.TEs_*_anchors.bed.gz
#boundary=2_dynamic/danRer10.TEs_*_boundaries.ctcf.bed.gz
#anchor=2_dynamic/danRer10.TEs_*_anchors.ctcf.bed.gz
#boundary=2_dynamic/danRer10.TEs_{brain,muscle}_boundaries.ctcf_atac.bed.gz
#anchor=2_dynamic/danRer10.TEs_{brain,muscle}_anchors.ctcf_atac.bed.gz
boundary=2_dynamic/danRer10.TEs_*_boundaries.atac_ctcf_zebrafish.bed.gz
anchor=2_dynamic/danRer10.TEs_*_anchors.cres_ctcf_zebrafish.bed.gz


# OUTPUT
dir_out=2_dynamic
out_cre=${dir_out}/danRer10.TEs_11_tissues_CREs.bed.gz
out_atac=${dir_out}/danRer10.TEs_11_tissues_ATAC.bed.gz
out_umr_lmr=${dir_out}/danRer10.TEs_11_tissues_UMR_LMR.bed.gz
out_dname=${dir_out}/danRer10.TE_wCG.DNAme_11_tissues_cov5.bed.gz
#out_boundary=${dir_out}/danRer10.TEs_2_tissues_boundaries.bed.gz
#out_anchor=${dir_out}/danRer10.TEs_2_tissues_anchors.bed.gz
#out_boundary=${dir_out}/danRer10.TEs_2_tissues_boundaries.ctcf.bed.gz
#out_anchor=${dir_out}/danRer10.TEs_2_tissues_anchors.ctcf.bed.gz
#out_boundary=${dir_out}/danRer10.TEs_2_tissues_boundaries.ctcf_atac.bed.gz
#out_anchor=${dir_out}/danRer10.TEs_2_tissues_anchors.ctcf_atac.bed.gz
out_boundary=${dir_out}/danRer10.TEs_2_tissues_boundaries.atac_ctcf_zebrafish.bed.gz
out_anchor=${dir_out}/danRer10.TEs_2_tissues_anchors.cres_ctcf_zebrafish.bed.gz


# COMMANDS
#gunzip $atac
paste $atac | awk -F"\t" -vOFS="\t" '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1,$2,$3,$4,$5,$6,$7); for (i=8;i<=NF;i+=8) {printf("\t%s", $i) } print ""}' | gzip -nc > $out_atac

#gunzip $cre
paste $cre | awk -F"\t" -vOFS="\t" '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1,$2,$3,$4,$5,$6,$7); for (i=8;i<=NF;i+=8) {printf("\t%s", $i) } print "" }' | gzip -nc > $out_cre

#gunzip $umr_lmr
paste $umr_lmr | awk -F"\t" -vOFS="\t" '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1,$2,$3,$4,$5,$6,$7); for (i=8;i<=NF;i+=8) {printf("\t%s", $i) } print ""}' | gzip -nc > $out_umr_lmr

gzip ${atac%.gz}
gzip ${cre%.gz}
gzip ${umr_lmr%.gz}

gunzip $dname
paste ${dname%.gz} | awk -F"\t" -vOFS="\t" '{printf("%s\t%s\t%s", $1,$2,$3); for (i=5;i<=NF;i+=7) {printf("\t%s", $i) } print "" }' | gzip -nc > $out_dname
gzip ${dname%.gz}

gunzip $boundary
paste ${boundary%.gz} | awk -F"\t" -vOFS="\t" '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1,$2,$3,$4,$5,$6,$7); for (i=8;i<=NF;i+=8) {printf("\t%s", $i) } print "" }' | gzip -nc > $out_boundary
gzip ${boundary%.gz}

gunzip $anchor
paste ${anchor%.gz} | awk -F"\t" -vOFS="\t" '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1,$2,$3,$4,$5,$6,$7); for (i=8;i<=NF;i+=8) {printf("\t%s", $i) } print "" }' | gzip -nc > $out_anchor
gzip ${anchor%.gz}




