#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID


# TE classes
bed_te_cg=danRer10.TE_wCG.bed.gz

# INPUT DATA
dir_in=0_cres
te_me=1_annot/danRer10.TE_wCG.DNAme_*_cov5.bed.gz


# OUTPUT
dir_out=2_dynamic
out=${dir_out}/te_dynamic_dname_across_10tissues.txt


# COMMANDS
intersectBed -sorted -wo -f 1 -F 1 -a $bed_te_cg -b $te_me |
    groupBy -c 13 -o collapse | sed "s/,/\t/g" |
    grep -w "NA" |
    awk -F"\t" 'BEGIN{ na=lo=me=hi=0 } 
        { for(i=4; i<=13; i++) {
            if( $i=="NA" ) { na++ }
            if( $i!="NA" && $i<0.25 ) { lo++ }
            if( $i!="NA" && $i>=0.25 && $i<0.75 ) { me++ }
            if( $i!="NA" && $i>=0.75 ) { hi++ } } }
        END{OFS="\t"; print "Missing", na, lo, me, hi } ' > $out

intersectBed -sorted -wo -f 1 -F 1 -a $bed_te_cg -b $te_me |
    groupBy -c 13 -o collapse | sed "s/,/\t/g" |
    awk -F"\t" 'BEGIN{ OFS="\t"} { for(i=4; i<=13; i++) {
            if( $i!="NA" && $i<0.25 ) { $i="LO" }
            else if( $i!="NA" && $i>=0.25 && $i<0.75 ) { $i="ME" }
            else if( $i!="NA" && $i>=0.75 ) { $i="HI" } } print $0}' |
    grep -w "LO" |
    awk -F"\t" 'BEGIN{ na=lo=me=hi=0 }
        { for(i=4; i<=13; i++) {
            if( $i=="NA" ) { na++ }
            if( $i=="LO" ) { lo++ }
            if( $i=="ME" ) { me++ }
            if( $i=="HI" ) { hi++ } } }
        END{OFS="\t"; print "Low", na, lo, me, hi } ' >> $out

intersectBed -sorted -wo -f 1 -F 1 -a $bed_te_cg -b $te_me |
    groupBy -c 13 -o collapse | sed "s/,/\t/g" |
    awk -F"\t" 'BEGIN{ OFS="\t"} { for(i=4; i<=13; i++) {
            if( $i!="NA" && $i<0.25 ) { $i="LO" }
            else if( $i!="NA" && $i>=0.25 && $i<0.75 ) { $i="ME" }
            else if( $i!="NA" && $i>=0.75 ) { $i="HI" } } print $0}' |
   grep -w "ME" |
    awk -F"\t" 'BEGIN{ na=lo=me=hi=0 }
        { for(i=4; i<=13; i++) {
            if( $i=="NA" ) { na++ }
            if( $i=="LO" ) { lo++ }
            if( $i=="ME" ) { me++ }
            if( $i=="HI" ) { hi++ } } }
        END{OFS="\t"; print "Intermediate", na, lo, me, hi } ' >> $out

intersectBed -sorted -wo -f 1 -F 1 -a $bed_te_cg -b $te_me |
    groupBy -c 13 -o collapse | sed "s/,/\t/g" |
    awk -F"\t" 'BEGIN{ OFS="\t"} { for(i=4; i<=13; i++) {
            if( $i!="NA" && $i<0.25 ) { $i="LO" }
            else if( $i!="NA" && $i>=0.25 && $i<0.75 ) { $i="ME" }
            else if( $i!="NA" && $i>=0.75 ) { $i="HI" } } print $0}' |
    grep -w "HI" |
    awk -F"\t" 'BEGIN{ na=lo=me=hi=0 }
        { for(i=4; i<=13; i++) {
            if( $i=="NA" ) { na++ }
            if( $i=="LO" ) { lo++ }
            if( $i=="ME" ) { me++ }
            if( $i=="HI" ) { hi++ } } }
        END{OFS="\t"; print "High", na, lo, me, hi } ' >> $out

