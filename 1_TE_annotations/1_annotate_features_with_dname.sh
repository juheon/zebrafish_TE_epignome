#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=5G
#SBATCH --array=1-13 #14
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID

# GENOME FEATURES
list=features.txt
feature=$( cat $list | sed "${ID}q;d" )
bed_feature=${feature}.bed.gz


# INPUT DATA
dss_files=fylab_WGBS_zt_*.dss.txt.gz

# OUTPUT
dir_out=1_annot
feature=${feature##*/}
feature=${feature#danRer10.}
feature=${feature#GRCz10.91.}
out=${dir_out}/${feature}.annotation_dname.txt

# COMMANDS
echo $feature > $out

## For each feature
lo=0
me=0
hi=0
na=0

for dss in $dss_files
do 
    zcat $dss | sed "/X$/d" |
    awk -vOFS="\t" '{print $1,$2,$2+2,$3,$4}' |
    intersectBed -sorted -f 1 -a stdin -b $bed_feature |
    awk 'BEGIN {lo=0; me=0; hi=0; na=0} $4<5 { na++ } 
        $4>=5 && $5/$4<0.25 { lo++ }
        $4>=5 && $5/$4>=0.25 && $5/$4<0.75 { me++ }
        $4>=5 && $5/$4>=0.75 { hi++ }
        END {print na"\n"lo"\n"me"\n"hi}' >> $out
done

