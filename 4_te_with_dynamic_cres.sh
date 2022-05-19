#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=dynamic_cre

# INPUT DATA
dir_in=0_cres
te_cre=2_dynamic/danRer10.TEs_11_tissues_CREs.bed.gz


# OUTPUT
dir_out=2_dynamic
out=${dir_out}/te_dynamic_cre_across_11_tissues.txt


# COMMANDS
zcat $te_cre | grep -w "." |
    awk -F"\t" 'BEGIN{ na_cnt=0; na=prom=weak=enhc=hete=0 } 
        { n=p=w=e=h=0;
	  for(i=8; i<=NF; i++) {
            if( $i=="." ) { n++; na_cnt++ }
            else if( $i=="Active_promoter" ) { p++ }
            else if( $i=="Weak_promoter" ) { w++ }
            else if( $i=="Active_enhancer" ) { e++ }
            else if( $i=="Heterochromatin" ) { h++ } }
	  na+=n*(n-1)/(n+p+w+e+h-1); prom+=n*p/(n+p+w+e+h-1); weak+=n*w/(n+p+w+e+h-1);
	  enhc+=n*e/(n+p+w+e+h-1); hete+=n*h/(n+p+w+e+h-1) }
        END{ printf ("No CRE\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\n", na, prom, weak, enhc, hete, na_cnt) } ' > $out

zcat $te_cre | grep -w "Active_promoter" |
    awk -F"\t" 'BEGIN{ prom_cnt=0; na=prom=weak=enhc=hete=0 } 
        { n=p=w=e=h=0;
	  for(i=8; i<=NF; i++) {
            if( $i=="." ) { n++ }
            else if( $i=="Active_promoter" ) { p++; prom_cnt++ }
            else if( $i=="Weak_promoter" ) { w++ }
            else if( $i=="Active_enhancer" ) { e++ }
            else if( $i=="Heterochromatin" ) { h++ } }
	  na+=p*n/(n+p+w+e+h-1); prom+=p*(p-1)/(n+p+w+e+h-1); weak+=p*w/(n+p+w+e+h-1);
	  enhc+=p*e/(n+p+w+e+h-1); hete+=p*h/(n+p+w+e+h-1) }
        END{ printf ("Active promoter\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\n", na, prom, weak, enhc, hete, prom_cnt) } ' >> $out

zcat $te_cre | grep -w "Weak_promoter" |
    awk -F"\t" 'BEGIN{ weak_cnt=0; na=prom=weak=enhc=hete=0 } 
        { n=p=w=e=h=0;
	  for(i=8; i<=NF; i++) {
            if( $i=="." ) { n++ }
            else if( $i=="Active_promoter" ) { p++ }
            else if( $i=="Weak_promoter" ) { w++; weak_cnt++ }
            else if( $i=="Active_enhancer" ) { e++ }
            else if( $i=="Heterochromatin" ) { h++ } }
	  na+=w*n/(n+p+w+e+h-1); prom+=w*p/(n+p+w+e+h-1); weak+=w*(w-1)/(n+p+w+e+h-1);
	  enhc+=w*e/(n+p+w+e+h-1); hete+=w*h/(n+p+w+e+h-1) }
        END{ printf ("Weak promoter\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\n", na, prom, weak, enhc, hete, weak_cnt) } ' >> $out

zcat $te_cre | grep -w "Active_enhancer" |
    awk -F"\t" 'BEGIN{ enhc_cnt=0; na=prom=weak=enhc=hete=0 } 
        { n=p=w=e=h=0;
	  for(i=8; i<=NF; i++) {
            if( $i=="." ) { n++ }
            else if( $i=="Active_promoter" ) { p++ }
            else if( $i=="Weak_promoter" ) { w++ }
            else if( $i=="Active_enhancer" ) { e++; enhc_cnt++ }
            else if( $i=="Heterochromatin" ) { h++ } }
	  na+=e*n/(n+p+w+e+h-1); prom+=e*p/(n+p+w+e+h-1); weak+=e*w/(n+p+w+e+h-1);
	  enhc+=e*(e-1)/(n+p+w+e+h-1); hete+=e*h/(n+p+w+e+h-1) }
        END{ printf ("Active enhancer\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\n", na, prom, weak, enhc, hete, enhc_cnt) } ' >> $out

zcat $te_cre | grep -w "Heterochromatin" |
    awk -F"\t" 'BEGIN{ hete_cnt=0; na=prom=weak=enhc=hete=0 } 
        { n=p=w=e=h=0;
	  for(i=8; i<=NF; i++) {
            if( $i=="." ) { n++ }
            else if( $i=="Active_promoter" ) { p++ }
            else if( $i=="Weak_promoter" ) { w++ }
            else if( $i=="Active_enhancer" ) { e++ }
            else if( $i=="Heterochromatin" ) { h++; hete_cnt++ } }
	  na+=h*n/(n+p+w+e+h-1); prom+=h*p/(n+p+w+e+h-1); weak+=h*w/(n+p+w+e+h-1);
	  enhc+=h*e/(n+p+w+e+h-1); hete+=h*(h-1)/(n+p+w+e+h-1) }
        END{ printf ("Heterochromatin\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\t%-10.2f\n", na, prom, weak, enhc, hete, hete_cnt) } ' >> $out

