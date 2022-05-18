#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=2,8
#SBATCH --job-name=hic


# SOFTWARE
module load bedtools/2.27.1

ID=$SLURM_ARRAY_TASK_ID
workdir=/scratch/twlab/hlee/zf_te


### INPUT
tissue_list=${workdir}/tissues-e.txt
dir_in=${workdir}/0_cres

tissue=$( cat $tissue_list | sed "${ID}q;d" )

tad=${dir_in}/TADs/TADs_${tissue}.bed
loop=${dir_in}/loop_anchors/loops_${tissue}.txt


### OUTPUT
boundaries=${dir_in}/TADs/TAD_boundaries_${tissue}.bed.gz
anchors=${dir_in}/loop_anchors/loop_anchors_${tissue}.bed.gz



### Commands
# TAD
cat $tad | sed "s/\r//g" |
    awk 'NR>1 && $1==chr && $2-1!=p {print chr"\t"p"\t"$2-1} {chr=$1; p=$3}' |
    sort -k1,1 -k2,2n | gzip -nc > $boundaries


# loops
cat $loop | awk '$0!~"^#" {print "chr"$1"\t"$2"\t"$3"\nchr"$4"\t"$5"\t"$6 }' |
    sort -u -k1,1 -k2,2n | gzip -nc > $anchors

