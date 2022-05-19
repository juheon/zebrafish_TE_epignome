#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=1-11
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID

# GENOME FEATURES
te=danRer10.TE.bed.gz
dna=danRer10.DNA.bed.gz
ltr=danRer10.LTR.bed.gz
line=danRer10.LINE.bed.gz
sine=danRer10.SINE.bed.gz
rc=danRer10.RC.bed.gz

te_bed=$(ls ${dir_te}/danRer10.{TE,DNA,LTR,LINE,SINE,RC}.bed.gz )


# INPUT DATA
list=tissues-e.txt
tissue=$( cat $list | sed "${ID}q;d" )

dss=fylab_WGBS_zt_${tissue}.dss.txt.gz

# OUTPUT
dir_out=1_annot
out=${dir_out}/TE_contribution_CGme_states_${tissue}.txt



# COMMANDS

for bed in $te_bed
do

    name=${bed##*/}
    name=${name%.bed.gz}
    name=${name#danRer10.}

    cnt_lo_te=$( zcat $dss | sed 1d | awk -v OFS="\t" '$3>=5 && $4/$3<0.25 {print $1,$2,$2+2}' | intersectBed -u -f 1 -sorted -a stdin -b $bed | wc -l )
    cnt_lo=$( zcat $dss | sed 1d | awk '$3>=5 && $4/$3<0.25' | wc -l )
    rate_lo=$( awk -v overlap="$cnt_lo_te" -v total="$cnt_lo" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_lo}\t${cnt_lo_te}\t${cnt_lo}\tLow\t${name}\t$tissue"

    cnt_me_te=$( zcat $dss | sed 1d | awk -v OFS="\t" '$3>=5 && $4/$3>=0.25 && $4/$3<0.75 {print $1,$2,$2+2}' | intersectBed -u -f 1 -sorted -a stdin -b $bed | wc -l )
    cnt_me=$( zcat $dss | sed 1d | awk '$3>=5 && $4/$3>=0.25 && $4/$3<0.75' | wc -l )
    rate_me=$( awk -v overlap="$cnt_me_te" -v total="$cnt_me" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_me}\t${cnt_me_te}\t${cnt_me}\tIntermediate\t${name}\t$tissue"

    cnt_hi_te=$( zcat $dss | sed 1d | awk -v OFS="\t" '$3>=5 && $4/$3>=0.75 {print $1,$2,$2+2}' | intersectBed -u -f 1 -sorted -a stdin -b $bed | wc -l )
    cnt_hi=$( zcat $dss | sed 1d | awk '$3>=5 && $4/$3>=0.75' | wc -l )
    rate_hi=$( awk -v overlap="$cnt_hi_te" -v total="$cnt_hi" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_hi}\t${cnt_hi_te}\t${cnt_hi}\tHigh\t${name}\t$tissue"

    cnt_na_te=$( zcat $dss | sed 1d | awk -v OFS="\t" '$3<5 {print $1,$2,$2+2}' | intersectBed -u -f 1 -sorted -a stdin -b $bed | wc -l )
    cnt_na=$( zcat $dss | sed 1d | awk '$3<5' | wc -l )
    rate_na=$( awk -v overlap="$cnt_na_te" -v total="$cnt_na" 'BEGIN{ printf("%.4f", overlap/total) }' )
    echo -e "${rate_na}\t${cnt_na_te}\t${cnt_na}\tMissing\t${name}\t$tissue"

done > $out


