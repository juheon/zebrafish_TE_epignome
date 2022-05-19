#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --array=1-11
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=cal_me_te
##SBATCH --mail-type=ALL

ID=$SLURM_ARRAY_TASK_ID

# SOFTWARES
module load bedtools/2.27.1

minCov=5
minCG=1


# INPUT
bed_te=/scratch/twlab/hlee/genomes/danRer10/rmsk/danRer10.TE_wCG.bed.gz

list=tissues-e.txt
tissue=$( cat $list | sed "${ID}q;d" )

dir_in=/scratch/twlab/hlee/fylab/wgbs/processing/4_dss
dss=${dir_in}/fylab_WGBS_zt_${tissue}.dss.txt.gz

# OUTPUT
dir_out=2_dynamic
out_me_te=${dir_out}/danRer10.TE_wCG.DNAme_${tissue}_cov5.bed.gz

# COMMANDS
zcat $dss | awk -vOFS="\t" 'NR>1 {print $1,$2,$2+2,$3,$4}' |
    intersectBed -sorted -F 1 -wo -a $bed_te -b stdin |
    awk -v min=$minCov -v OFS="\t" '$11>=min {print $1,$2,$3,$4,sprintf("%.4g",$12/$11),$6,$7}' |
    groupBy -c 5,5 -o mean,count -prec 4 -full |
    awk -v min=$minCG -v OFS="\t" '$9>=min {print $1,$2,$3,$4,$8,$6,$7}' |
    intersectBed -sorted -f 1 -F 1 -wao -a $bed_te -b stdin |
    awk -v OFS="\t" '$15==0 {$12="NA"} {print $1,$2,$3,$4,$12,$6,$7}' |
    gzip -nc  > $out_me_te

