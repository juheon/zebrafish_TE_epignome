#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
#ml htslib/1.3.1
#ml python2/2.7
ml R/3.6.1

aggregate_annot=../TEprof_zebrafish/bin/aggregateProcessedAnnotation.R
inputdir=../output_teprof_zebrafish/03_annotation_tpm_process_output
outputdir=../output_teprof_zebrafish/04_aggregate_annotation_process_output
argument=`pwd`/arguments.txt

mkdir $outputdir

### INPUT # 26 annotated all_c
#c_files=$( ls ${workdir}/*_rep1.gtf_annotated_filtered_test_all_c )

log=$inputdir/aggregateProcessedAnnotation.R.out

$aggregate_annot -a $argument -i $inputdir &>$log

wait 
mv $log $outputdir
mv $inputdir/initial_candidate_list.tsv $outputdir
mv $inputdir/filter_combined_candidates.tsv $outputdir
mv $inputdir/Step4.RData $outputdir


