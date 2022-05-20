#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
ml htslib/1.3.1
ml python2/2.7
ml cufflinks/2.2.1

command_update=`pwd`/../TEprof_zebrafish/bin/rmskhg38_annotate_gtf_update_test_tpm_cuff.py

argument=`pwd`/arguments.txt
inputdir=`pwd`/../output_teprof_zebrafish/07_merge_ref_gtf
outputdir=`pwd`/../output_teprof_zebrafish/08_annot_merged_gtf

### INPUT


### COMMANDS
$command_update $inputdir/reference_merged_candidates.gff3 $argument

mv $inputdir/*_test_all $outputdir



