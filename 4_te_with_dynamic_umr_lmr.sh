#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=dynamic_umr_lmr

# INPUT DATA
dir_in=0_cres
te_umr_lmr=2_dynamic/danRer10.TEs_11_tissues_UMR_LMR.bed.gz


# OUTPUT
dir_out=2_dynamic
out=${dir_out}/te_dynamic_umr_lmr_across_11_tissues.txt


# COMMANDS
zcat $te_umr_lmr | grep -w "." |
    awk -F"\t" 'BEGIN{ na_cnt=0; na=umr=lmr=0 } 
        { n=u=l=0;
	  for(i=8; i<=NF; i++) {
            if( $i=="." ) { n++; na_cnt++ }
            else if( $i=="UMR" ) { u++ }
            else if( $i=="LMR" ) { l++ } }
	  na+=n*(n-1)/(n+u+l-1); umr+=n*u/(n+u+l-1); lmr+=n*l/(n+u+l-1) }
        END{ printf("No CRE\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\n", na, umr, lmr, na_cnt) } ' > $out

zcat $te_umr_lmr | grep -w "UMR" |
    awk -F"\t" 'BEGIN{ umr_cnt=0; na=umr=lmr=0 } 
        { n=u=l=0;
	  for(i=8; i<=NF; i++) {
            if( $i=="." ) { n++ }
            else if( $i=="UMR" ) { u++; umr_cnt++ }
            else if( $i=="LMR" ) { l++ } }
	  na+=u*n/(n+u+l-1); umr+=u*(u-1)/(n+u+l-1); lmr+=u*l/(n+u+l-1) }
        END{ printf("UMR\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\n", na, umr, lmr, umr_cnt) } ' >> $out

zcat $te_umr_lmr | grep -w "LMR" |
    awk -F"\t" 'BEGIN{ lmr_cnt=0; na=umr=lmr=0 } 
        { n=u=l=0;
	  for(i=8; i<=NF; i++) {
            if( $i=="." ) { n++ }
            else if( $i=="UMR" ) { u++ }
            else if( $i=="LMR" ) { l++; lmr_cnt++ } }
	  na+=l*n/(n+u+l-1); umr+=l*u/(n+u+l-1); lmr+=l*(l-1)/(n+u+l-1) }
        END{ printf("LMR\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\n", na, umr, lmr, lmr_cnt) } ' >> $out

