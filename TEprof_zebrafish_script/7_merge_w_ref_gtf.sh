#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
ml htslib/1.3.1
ml python2/2.7
ml cufflinks/2.2.1

argument=./arguments.txt

ref=`pwd`/../annotation/Danio_rerio.GRCz10.91.gtf
inputdir=../output_teprof_zebrafish/06_filter_candi_read_info
outputdir=../output_teprof_zebrafish/07_merge_ref_gtf

mkdir $outputdir
### INPUT


### COMMANDS
#gffread -E $inputdir/candidate_transcripts.gff3 -T -o $outputdir/candidate_transcripts.gtf

cd $outputdir
#echo candidate_transcripts.gtf > cuffmergegtf.list


# cufflinks version 2.2.1 (cuffmerge = merge_cuff_asms v 1.0.0)
# taco might be better

log=cuffmerge.log

cuffmerge -o ./merged_asm_full -g $ref cuffmergegtf.list &> $log
wait
mv ./merged_asm_full/merged.gtf reference_merged_candidates.gtf

gffread -E reference_merged_candidates.gtf -o- > reference_merged_candidates.gff3
