#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
##SBATCH --array=1-10
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=share_tissue

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID


# TE classes
dir_te=/scratch/twlab/hlee/genomes/danRer10/rmsk
bed_te_cg=${dir_te}/danRer10.TE_wCG.bed.gz

# INPUT DATA
te_me=2_dynamic/danRer10.TE_wCG.DNAme_*_cov5.bed.gz


# OUTPUT
dir_out=2_dynamic
out_te_me=${dir_out}/cnt_te_me_num_shared_tissue.txt


# COMMANDS
intersectBed -sorted -wo -f 1 -F 1 -a $bed_te_cg -b $te_me |
    groupBy -c 13 -o collapse | sed "s/,/\t/g" |
    awk -F"\t" 'BEGIN{ for (i=0; i<=10; i++) {cnt_na[i]=0; cnt_lo[i]=0; cnt_me[i]=0; cnt_hi[i]=0 } }
        { na=0; lo=0; me=0; hi=0;
        for(i=4; i<=NF; i++) {
            if( $i=="NA" ) {na++}
            if( $i!="NA" && $i<0.25 ) {lo++}
            if( $i!="NA" && $i>=0.25 && $i<0.75 ) {me++}
            if( $i!="NA" && $i>=0.75 ) {hi++} }
        cnt_na[na]++; cnt_lo[lo]++; cnt_me[me]++; cnt_hi[hi]++ }
        END{ OFS="\t"; print "NA", "Low", "Intermediate", "High";
        for(i=1; i<=11; i++) { print cnt_na[i],cnt_lo[i],cnt_me[i],cnt_hi[i] } } ' > $out_te_me

