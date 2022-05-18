#!/bin/bash
# Author: Hyung Joo Lee


# SOFTWARE
ml bedtools/2.27.1

# GENOME FEATURES
subfam_bed=$(ls ../annotation/TEsubfamily/*bed.gz)
genome_size=../annotation/danRer10.chrom.sizes


# INPUT DATA
list=../tissues-e.txt
tissue_list=($( cat $list  ))

for i in ${tissue_list[@]}
do
	tissue=$i
	bed_prom=../0_CRE_file_prepared/tissue_specific/${tissue}_active_promoter.tissue-specific.bed.gz
	bed_weak=../0_CRE_file_prepared/tissue_specific/${tissue}_weak_promoter.tissue-specific.bed.gz
	bed_enhc=../0_CRE_file_prepared/tissue_specific/${tissue}_active_enhancer.tissue-specific.bed.gz
	
	bed_prox=../0_CRE_file_prepared/tissue_specific/${tissue}_proximal_ATAC.tissue-specific.bed.gz
	bed_dist=../0_CRE_file_prepared/tissue_specific/${tissue}_distal_ATAC.tissue-specific.bed.gz
	
	bed_umr=../0_CRE_file_prepared/tissue_specific/${tissue}_UMR.tissue-specific.bed.gz
	bed_lmr=../0_CRE_file_prepared/tissue_specific/${tissue}_LMR.tissue-specific.bed.gz
	
	
	# OUTPUT
	dir_out=`pwd`	
	out_prom=${dir_out}/subfam_enriched_active_promoter_${tissue}.txt
	out_weak=${dir_out}/subfam_enriched_weak_promoter_${tissue}.txt
	out_enhc=${dir_out}/subfam_enriched_enhancer_${tissue}.txt
	
	out_prox=${dir_out}/subfam_enriched_proximal_ATAC_${tissue}.txt
	out_dist=${dir_out}/subfam_enriched_distal_ATAC_${tissue}.txt
	
	out_umr=${dir_out}/subfam_enriched_UMR_${tissue}.txt
	out_lmr=${dir_out}/subfam_enriched_LMR_${tissue}.txt
	
	
	# COMMANDS
	len_genome=$( awk 'BEGIN{s=0} {s+=$2} END{print s}' $genome_size )
	
	len_prom=$( mergeBed -i $bed_prom | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	len_weak=$( mergeBed -i $bed_weak | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	len_enhc=$( mergeBed -i $bed_enhc | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	
	len_prox=$( mergeBed -i $bed_prox | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	len_dist=$( mergeBed -i $bed_dist | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	
	len_umr=$( mergeBed -i $bed_umr | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	len_lmr=$( mergeBed -i $bed_lmr | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	for bed_te in ${subfam_bed[@]}
	do
	    name_te=${bed_te##*/}
	    name_te=${name_te%.bed.gz}
	    len_te=$( mergeBed -i $bed_te | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	
	    bp_prom_te=$( intersectBed -sorted -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_prom | sort -k1,1 -k2,2n | mergeBed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	    lor=$( awk -v nn=$bp_prom_te -v nd=$len_te -v dn=$len_prom -v dd=$len_genome 'BEGIN{ if (nn>0) { printf ("%.4f", log( (nn/nd)/(dn/dd) ) / log(2) ) } else {printf("NA")} }' )
	    cnt_prom_te=$( intersectBed -sorted -wa -u -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_prom | wc -l )
	    if [ $cnt_prom_te -le 10 ]; then
	        lor="NA"; fi
	    echo -e "$lor\t$name_te\t$tissue\tActive_promoter\t$bp_prom_te\t$len_te\t$len_prom" >> $out_prom
	  
	
	    bp_weak_te=$( intersectBed -sorted -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_weak | sort -k1,1 -k2,2n | mergeBed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	    lor=$( awk -v nn=$bp_weak_te -v nd=$len_te -v dn=$len_weak -v dd=$len_genome 'BEGIN{ if (nn>0) { printf ("%.4f", log( (nn/nd)/(dn/dd) ) / log(2) ) } else {printf("NA")} }' )
	    cnt_weak_te=$( intersectBed -sorted -wa -u -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_weak | wc -l )
	    if [ $cnt_weak_te -le 10 ]; then
	        lor="NA"; fi
	    echo -e "$lor\t$name_te\t$tissue\tWeak_promoter\t$bp_weak_te\t$len_te\t$len_weak" >> $out_weak
	
	    bp_enhc_te=$( intersectBed -sorted -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_enhc | sort -k1,1 -k2,2n | mergeBed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	    lor=$( awk -v nn=$bp_enhc_te -v nd=$len_te -v dn=$len_enhc -v dd=$len_genome 'BEGIN{ if (nn>0) { printf ("%.4f", log( (nn/nd)/(dn/dd) ) / log(2) ) } else {printf("NA")} }' )
	    cnt_enhc_te=$( intersectBed -sorted -wa -u -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_enhc | wc -l )
	    if [ $cnt_enhc_te -le 10 ]; then
	        lor="NA"; fi
	    echo -e "$lor\t$name_te\t$tissue\tActive_enhancer\t$bp_enhc_te\t$len_te\t$len_enhc" >> $out_enhc
	
	
	
	    bp_prox_te=$( intersectBed -sorted -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_prox | sort -k1,1 -k2,2n | mergeBed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	    lor=$( awk -v nn=$bp_prox_te -v nd=$len_te -v dn=$len_prox -v dd=$len_genome 'BEGIN{ if (nn>0) { printf ("%.4f", log( (nn/nd)/(dn/dd) ) / log(2) ) } else {printf("NA")} }' )
	    cnt_prox_te=$( intersectBed -sorted -wa -u -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_prox | wc -l )
	    if [ $cnt_prox_te -le 10 ]; then
	        lor="NA"; fi
	    echo -e "$lor\t$name_te\t$tissue\tProximal_ATAC\t$bp_prox_te\t$len_te\t$len_prox" >> $out_prox
	
	    bp_dist_te=$( intersectBed -sorted -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_dist | sort -k1,1 -k2,2n | mergeBed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	    lor=$( awk -v nn=$bp_dist_te -v nd=$len_te -v dn=$len_dist -v dd=$len_genome 'BEGIN{ if (nn>0) { printf ("%.4f", log( (nn/nd)/(dn/dd) ) / log(2) ) } else {printf("NA")} }' )
	    cnt_dist_te=$( intersectBed -sorted -wa -u -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_dist | wc -l )
	    if [ $cnt_dist_te -le 10 ]; then
	        lor="NA"; fi
	    echo -e "$lor\t$name_te\t$tissue\tDistal_ATAC\t$bp_dist_te\t$len_te\t$len_dist" >> $out_dist
	
	
	    bp_umr_te=$( intersectBed -sorted -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_umr | sort -k1,1 -k2,2n | mergeBed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	    lor=$( awk -v nn=$bp_umr_te -v nd=$len_te -v dn=$len_umr -v dd=$len_genome 'BEGIN{ if (nn>0) { printf ("%.4f", log( (nn/nd)/(dn/dd) ) / log(2) ) } else {printf("NA")} }' )
	    cnt_umr_te=$( intersectBed -sorted -wa -u -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_umr | wc -l )
	    if [ $cnt_umr_te -le 10 ]; then
	        lor="NA"; fi
	    echo -e "$lor\t$name_te\t$tissue\tUMR\t$bp_umr_te\t$len_te\t$len_umr" >> $out_umr
	
	
	    bp_lmr_te=$( intersectBed -sorted -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_lmr | sort -k1,1 -k2,2n | mergeBed | awk 'BEGIN{s=0} {s+=$3-$2} END{print s}' )
	    lor=$( awk -v nn=$bp_lmr_te -v nd=$len_te -v dn=$len_lmr -v dd=$len_genome 'BEGIN{ if (nn>0) { printf ("%.4f", log( (nn/nd)/(dn/dd) ) / log(2) ) } else {printf("NA")} }' )
	    cnt_lmr_te=$( intersectBed -sorted -wa -u -a <(zcat $bed_te | sort -k1,1 -k2,2n ) -b $bed_lmr | wc -l )
	    if [ $cnt_lmr_te -le 10 ]; then
	        lor="NA"; fi
	    echo -e "$lor\t$name_te\t$tissue\tLMR\t$bp_lmr_te\t$len_te\t$len_lmr" >> $out_lmr
		
	done
done

