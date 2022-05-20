#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
ml htslib/1.3.1
ml python2/2.7
ml stringtie/2.1.4

inputdir=`pwd`/../output_teprof_zebrafish/07_merge_ref_gtf
outputdir=`pwd`/../output_teprof_zebrafish/09_cal_trx_level_tpm
## This is the path in HTCF cluster
bam_dir=`pwd`/../stringtie_bam/

argument=./arguments.txt


### INPUT
tsv=`pwd`/.././04_aggregate_annotation_process_output/filter_combined_candidates.tsv
list=`pwd`/.././04_aggregate_annotation_process_output/samples.txt


### COMMANDS
flist=($(find $bam_dir -maxdepth 1 -mindepth 1 -name "*bam"))
for file in ${flist[@]}
do 
    xbase=${file##*/}
	echo $xbase
    echo -e "samtools view -q 255 -h $file | stringtie - -o $outputdir/${xbase%.*}.gtf -e -b $outputdir/${xbase%.*}_stats -p 2 -m 100 -c 1 -G $inputdir/reference_merged_candidates.gtf" >> $outputdir/quantificationCommands.txt
done

while read -r line
do
        eval $line &
done < $outputdir/quantificationCommands.txt

