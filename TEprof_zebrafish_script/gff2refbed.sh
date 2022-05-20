#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
ml kentUCSC/334
ml bedops/2.4.19
ml htslib/1.3.1

workdir=[]


### INPUT
in_gtf=${workdir}/07_merge_ref_gtf/reference_merged_candidates.gtf

### OUTPUT
out_genepred=${workdir}/reference_merged_candidates.genePred
out_refbed=${workdir}/reference_merged_candidates.refBed.gz

### COMMANDS
gtfToGenePred -genePredExt -geneNameAsName2 $in_gtf $out_genepred

cat $out_genepred | awk -F"\t" -vOFS="\t" '{print $2,$4,$5,$6,$7,$3,$12,$1,"coding",$9,$10,"."}' |
    sort -k1,1 -k2,2n | bgzip > $out_refbed
tabix -p bed $out_refbed


