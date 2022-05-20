#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
ml htslib/1.3.1

annot_gtf=../TEprof_zebrafish/bin/rmskhg38_annotate_gtf_update_test_tpm.py
inputdir=../output_teprof_zebrafish/01_input
outputdir=../output_teprof_zebrafish/02_annotate_each_gtf_output
argument=./arguments.txt

mkdir $outputdir

### INPUT # 26 samples 
gtf_files=$( ls $inputdir/*.gtf )


for gtf in $gtf_files
do

    ### OUTPUT
    log=${gtf/gtf/teprof_annot_gtf.log}

    ### COMMANDS
    mv $gtf $gtf.bak
    sed "/^chrUn/d" $gtf.bak > $gtf
    $annot_gtf $gtf $argument &>$log &
done

wait
mv $inputdir/*.teprof_annot_gtf.log $outputdir
mv $inputdir/*_all $outputdir 


