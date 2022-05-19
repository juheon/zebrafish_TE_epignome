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
tad=2_dynamic/danRer10.TEs_2_tissues_boundaries.bed.gz
loop=2_dynamic/danRer10.TEs_2_tissues_anchors.bed.gz

# OUTPUT
dir_out=2_dynamic
out_tad=${dir_out}/cnt_te_tad_shared_tissue.txt
out_loop=${dir_out}/cnt_te_loop_shared_tissue.txt

# COMMANDS
zcat $tad |
    awk -F"\t" 'BEGIN{ for (i=0; i<=2; i++) {cnt_tad[i]=0 } }
        { tad=0;
        for(i=4; i<=NF; i++) {
            if( $i=="TAD_boundary" ) {tad++} }
        cnt_tad[tad]++ }
        END{ OFS="\t"; print "TAD_boundary";
        for(i=0; i<=2; i++) { print cnt_tad[i] } } ' > $out_tad

zcat $loop |
    awk -F"\t" 'BEGIN{ for (i=0; i<=2; i++) {cnt_loop[i]=0 } }
        { loop=0;
        for(i=4; i<=NF; i++) {
            if( $i=="loop_anchor" ) {loop++} }
        cnt_loop[loop]++ }
        END{ OFS="\t"; print "Loop_anchor";
        for(i=0; i<=2; i++) { print cnt_loop[i] } } ' > $out_loop




