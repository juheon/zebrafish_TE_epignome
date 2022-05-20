#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
#ml htslib/1.3.1
#ml python2/2.7

ml R/3.6.1

r_script=`pwd`/../TEprof_zebrafish/bin/finalStatisticsOutput.R

argument=`pwd`/./arguments.txt
inputdir=`pwd`/../output_teprof_zebrafish/10_process_output
outputdir=`pwd`/../output_teprof_zebrafish/11_final_table


### INPUT # 26 annotated all_c
#c_files=$( ls ${workdir}/*_rep1.gtf_annotated_filtered_test_all_c )

log=$outputdir/finalStatisticsOutput.R.out

$r_script -a $argument -s $inputdir/Step10.RData -p $inputdir -o $outputdir &>$log

wait
mv $inputdir/all*tsv $outputdir 
mv $inputdir/All*xlsx $outputdir
mv $inputdir/merged*refBed $outputdir
mv $inputdir/Step11*RData $outputdir
