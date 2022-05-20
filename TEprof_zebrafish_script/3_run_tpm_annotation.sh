#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
ml htslib/1.3.1
ml python2/2.7

tpm_annot=../TEprof_zebrafish/bin/annotationtpmprocess.py


inputdir=../output_teprof_zebrafish/02_annotate_each_gtf_output
outputdir=../output_teprof_zebrafish/03_annotation_tpm_process_output
argument=./arguments.txt

mkdir $outputdir
### INPUT 26 files
annot_files=$( ls $inputdir/*.gtf_annotated_filtered_test_all )


for annot in $annot_files
do

    ### OUTPUT
    log=${annot/.gtf_annotated_filtered_test_all/.teprof_tpm_annot.log}

    ### COMMANDS
    python $tpm_annot $annot  &>$log &

done

wait 
mv $inputdir/*_tpm_annot.log $outputdir
mv $inputdir/*_all_c $outputdir


