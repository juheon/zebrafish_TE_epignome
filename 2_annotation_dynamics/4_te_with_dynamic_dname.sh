#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --job-name=dynamic_cre

# INPUT DATA
dir_in=0_cres
te_dname=2_dynamic/danRer10.TE_wCG.DNAme_11_tissues_cov5.bed.gz


# OUTPUT
dir_out=2_dynamic
out=${dir_out}/te_dynamic_dname_across_11_tissues.txt


# COMMANDS
zcat $te_dname |
    awk -F"\t" 'BEGIN{ OFS="\t"} { for(i=4; i<=NF; i++) {
            if( $i!="NA" && $i<0.25 ) { $i="LO" }
            else if( $i!="NA" && $i>=0.25 && $i<0.75 ) { $i="ME" }
            else if( $i!="NA" && $i>=0.75 ) { $i="HI" } } print $0}' |
    grep -w "NA" |
    awk -F"\t" 'BEGIN{ na_cnt=0; na=lo=me=hi=0 } 
        { n=l=m=h=0;
	  for(i=4; i<=NF; i++) {
            if( $i=="NA" ) { n++; na_cnt++ }
            else if( $i=="LO" ) { l++ }
            else if( $i=="ME" ) { m++ }
            else if( $i=="HI" ) { h++ } }
	  na+=n*(n-1)/(n+l+m+h-1); lo+=n*l/(n+l+m+h-1); me+=n*m/(n+l+m+h-1); hi+=n*h/(n+l+m+h-1) }
        END{ printf ("Missing\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\n", na, lo, me, hi, na_cnt) } ' > $out

zcat $te_dname |
    awk -F"\t" 'BEGIN{ OFS="\t"} { for(i=4; i<=NF; i++) {
            if( $i!="NA" && $i<0.25 ) { $i="LO" }
            else if( $i!="NA" && $i>=0.25 && $i<0.75 ) { $i="ME" }
            else if( $i!="NA" && $i>=0.75 ) { $i="HI" } } print $0}' |
    grep -w "LO" |
    awk -F"\t" 'BEGIN{ lo_cnt=0; na=lo=me=hi=0 } 
        { n=l=m=h=0;
	  for(i=4; i<=NF; i++) {
            if( $i=="NA" ) { n++ }
            else if( $i=="LO" ) { l++; lo_cnt++ }
            else if( $i=="ME" ) { m++ }
            else if( $i=="HI" ) { h++ } }
	  na+=l*n/(n+l+m+h-1); lo+=l*(l-1)/(n+l+m+h-1); me+=l*m/(n+l+m+h-1); hi+=l*h/(n+l+m+h-1) }
        END{ printf ("Low\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\n", na, lo, me, hi, lo_cnt) } ' >> $out

zcat $te_dname |
    awk -F"\t" 'BEGIN{ OFS="\t"} { for(i=4; i<=NF; i++) {
            if( $i!="NA" && $i<0.25 ) { $i="LO" }
            else if( $i!="NA" && $i>=0.25 && $i<0.75 ) { $i="ME" }
            else if( $i!="NA" && $i>=0.75 ) { $i="HI" } } print $0}' |
    grep -w "ME" |
    awk -F"\t" 'BEGIN{ me_cnt=0; na=lo=me=hi=0 } 
        { n=l=m=h=0;
	  for(i=4; i<=NF; i++) {
            if( $i=="NA" ) { n++ }
            else if( $i=="LO" ) { l++ }
            else if( $i=="ME" ) { m++; me_cnt++ }
            else if( $i=="HI" ) { h++ } }
	  na+=m*n/(n+l+m+h-1); lo+=m*l/(n+l+m+h-1); me+=m*(m-1)/(n+l+m+h-1); hi+=m*h/(n+l+m+h-1) }
        END{ printf ("Intermediate\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\n", na, lo, me, hi, me_cnt) } ' >> $out

zcat $te_dname |
    awk -F"\t" 'BEGIN{ OFS="\t"} { for(i=4; i<=NF; i++) {
            if( $i!="NA" && $i<0.25 ) { $i="LO" }
            else if( $i!="NA" && $i>=0.25 && $i<0.75 ) { $i="ME" }
            else if( $i!="NA" && $i>=0.75 ) { $i="HI" } } print $0}' |
    grep -w "HI" |
    awk -F"\t" 'BEGIN{ hi_cnt=0; na=lo=me=hi=0 } 
        { n=l=m=h=0;
	  for(i=4; i<=NF; i++) {
            if( $i=="NA" ) { n++ }
            else if( $i=="LO" ) { l++ }
            else if( $i=="ME" ) { m++ }
            else if( $i=="HI" ) { h++; hi_cnt++ } }
	  na+=h*n/(n+l+m+h-1); lo+=h*l/(n+l+m+h-1); me+=h*m/(n+l+m+h-1); hi+=h*(h-1)/(n+l+m+h-1) }
        END{ printf ("High\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\n", na, lo, me, hi, hi_cnt) } ' >> $out


