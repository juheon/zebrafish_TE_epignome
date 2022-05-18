#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=1-11
#SBATCH --job-name=queiscnt


# SOFTWARE
module load bedtools/2.27.1

ID=$SLURM_ARRAY_TASK_ID
workdir=/scratch/twlab/hlee/zf_te


### INPUT
chrSize=/scratch/twlab/hlee/genomes/danRer10/danRer10.chrom.sizes

tissue_list=${workdir}/tissues-e.txt
tissue=$( cat $tissue_list | sed "${ID}q;d" )

dir_in=${workdir}/0_cres
prom=${dir_in}/activepromoter/${tissue}.active.promoter
weak=${dir_in}/weakpromoter/${tissue}.weak.promoter
enhc=${dir_in}/enhancer/${tissue}.active.enhancer
hete=${dir_in}/heterochromatin/${tissue}.Hetero-chromatin.bed


### OUTPUT
dir_out=${workdir}/0_cres/queiscent
queiscent=${dir_out}/queiscent_${tissue}.bed.gz


### Commands
sort -k1,1 $chrSize |
    awk -vOFS="\t" '{print $1,0,$2}' |
    subtractBed -a stdin -b $prom $weak $enhc $hete |
    gzip -nc > $queiscent

