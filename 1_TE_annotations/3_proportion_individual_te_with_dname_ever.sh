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
bed_te_cg=${dir_te}/danRer10.TE_wCG.bed.gz

# INPUT DATA
te_me=2_dynamic/danRer10.TE_wCG.DNAme_*_cov5.bed.gz


# OUTPUT
dir_out=2_dynamic
out1=${dir_out}/fraction_dname_states_ever_te_class.txt
out2=${dir_out}/fraction_dname_states_ever.txt

# COMMANDS
cnt_te_cg=$( zcat $bed_te_cg | wc -l )

na_ever=$( zcat $te_me | awk '$5=="NA"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
na_dna=$( zcat $te_me | awk '$5=="NA" && $7=="DNA"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
na_ltr=$( zcat $te_me | awk '$5=="NA" && $7=="LTR"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
na_line=$( zcat $te_me | awk '$5=="NA" && $7=="LINE"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
na_sine=$( zcat $te_me | awk '$5=="NA" && $7=="SINE"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
na_rc=$( zcat $te_me | awk '$5=="NA" && $7=="RC"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )

lo_ever=$( zcat $te_me | awk '$5!="NA" && $5<0.25' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
lo_dna=$( zcat $te_me | awk '$5!="NA" && $5<0.25 && $7=="DNA"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
lo_ltr=$( zcat $te_me | awk '$5!="NA" && $5<0.25 && $7=="LTR"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
lo_line=$( zcat $te_me | awk '$5!="NA" && $5<0.25 && $7=="LINE"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
lo_sine=$( zcat $te_me | awk '$5!="NA" && $5<0.25 && $7=="SINE"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
lo_rc=$( zcat $te_me | awk '$5!="NA" && $5<0.25 && $7=="RC"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )

me_ever=$( zcat $te_me | awk '$5!="NA" && $5>=0.25 && $5<0.75' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
me_dna=$( zcat $te_me | awk '$5!="NA" && $5>=0.25 && $5<0.75 && $7=="DNA"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
me_ltr=$( zcat $te_me | awk '$5!="NA" && $5>=0.25 && $5<0.75 && $7=="LTR"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
me_line=$( zcat $te_me | awk '$5!="NA" && $5>=0.25 && $5<0.75 && $7=="LINE"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
me_sine=$( zcat $te_me | awk '$5!="NA" && $5>=0.25 && $5<0.75 && $7=="SINE"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
me_rc=$( zcat $te_me | awk '$5!="NA" && $5>=0.25 && $5<0.75 && $7=="RC"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )

hi_ever=$( zcat $te_me | awk '$5!="NA" && $5>=0.75' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
hi_dna=$( zcat $te_me | awk '$5!="NA" && $5>=0.75 && $7=="DNA"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
hi_ltr=$( zcat $te_me | awk '$5!="NA" && $5>=0.75 && $7=="LTR"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
hi_line=$( zcat $te_me | awk '$5!="NA" && $5>=0.75 && $7=="LINE"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
hi_sine=$( zcat $te_me | awk '$5!="NA" && $5>=0.75 && $7=="SINE"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )
hi_rc=$( zcat $te_me | awk '$5!="NA" && $5>=0.75 && $7=="RC"' | cut -f1-4 | sort -u -k1,1 -k2,2n | wc -l )

echo -e "Missing\t$na_dna\t$na_ltr\t$na_line\t$na_sine\t$na_rc" >$out1
echo -e "Low\t$lo_dna\t$lo_ltr\t$lo_line\t$lo_sine\t$lo_rc" >>$out1
echo -e "Medium\t$me_dna\t$me_ltr\t$me_line\t$me_sine\t$me_rc" >>$out1
echo -e "High\t$hi_dna\t$hi_ltr\t$hi_line\t$hi_sine\t$hi_rc" >>$out1

awk -v"a=$na_ever" -v"b=$cnt_te_cg" 'BEGIN{printf("%.5g\t%d\t%d\tMissing\n", a/b,a,b)}' >$out2
awk -v"a=$lo_ever" -v"b=$cnt_te_cg" 'BEGIN{printf("%.5g\t%d\t%d\tLow\n", a/b,a,b)}' >>$out2
awk -v"a=$me_ever" -v"b=$cnt_te_cg" 'BEGIN{printf("%.5g\t%d\t%d\tMedium\n", a/b,a,b)}' >>$out2
awk -v"a=$hi_ever" -v"b=$cnt_te_cg" 'BEGIN{printf("%.5g\t%d\t%d\tHigh\n", a/b,a,b)}' >>$out2

