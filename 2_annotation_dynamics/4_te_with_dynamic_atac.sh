#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --job-name=dynamic_atac

# INPUT DATA
dir_in=0_cres
te_atac=2_dynamic/danRer10.TEs_11_tissues_ATAC.bed.gz


# OUTPUT
dir_out=2_dynamic
out=${dir_out}/te_dynamic_atac_across_11_tissues.txt


# COMMANDS
zcat $te_atac | grep -w "." |
    awk -F"\t" 'BEGIN{ na_cnt=0; na=prox=dist=0 } 
        { n=p=d=0;
	  for(i=8; i<=NF; i++) {
            if( $i=="." ) { n++; na_cnt++ }
            else if( $i=="Proximal" ) { p++ }
            else if( $i=="Distal" ) { d++ } }
	  na+=n*(n-1)/(n+p+d-1); prox+=n*p/(n+p+d-1); dist+=n*d/(n+p+d-1) }
        END{ printf("No CRE\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\n", na, prox, dist, na_cnt) } ' > $out

zcat $te_atac | grep -w "Proximal" |
    awk -F"\t" 'BEGIN{ prox_cnt=0; na=prox=dist=0 } 
        { n=p=d=0;
	  for(i=8; i<=NF; i++) {
            if( $i=="." ) { n++ }
            else if( $i=="Proximal" ) { p++; prox_cnt++ }
            else if( $i=="Distal" ) { d++ } }
	  na+=p*n/(n+p+d-1); prox+=p*(p-1)/(n+p+d-1); dist+=p*d/(n+p+d-1) }
        END{ printf("Proximal\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\n", na, prox, dist, prox_cnt) } ' >> $out

zcat $te_atac | grep -w "Distal" |
    awk -F"\t" 'BEGIN{ dist_cnt=0; na=prox=dist=0 } 
        { n=p=d=0;
	  for(i=8; i<=NF; i++) {
            if( $i=="." ) { n++ }
            else if( $i=="Proximal" ) { p++ }
            else if( $i=="Distal" ) { d++; dist_cnt++ } }
	  na+=d*n/(n+p+d-1); prox+=d*p/(n+p+d-1); dist+=d*(d-1)/(n+p+d-1) }
        END{ printf("Distal\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\n", na, prox, dist, dist_cnt) } ' >> $out

