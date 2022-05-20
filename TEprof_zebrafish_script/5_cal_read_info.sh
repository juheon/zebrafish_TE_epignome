#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
ml htslib/1.3.1
ml python2/2.7

commandsmax_pe=`pwd`/../TEprof_zebrafish/bin/commandsmax_speed.py
commandsmax_se=`pwd`/../TEprof_zebrafish/bin/commandsmax_speed_se.py

inputdir=`pwd`/../output_teprof_zebrafish/04_aggregate_annotation_process_output
outputdir=`pwd`/../output_teprof_zebrafish/05_cal_read_info

bam_dir=`pwd`/../stringtie_bam/

### INPUT
tsv=$inputdir/filter_combined_candidates.tsv
list=$inputdir/samples.txt


### COMMANDS
mkdir $outputdir/filterreadstats

cd $inputdir
$commandsmax_pe $tsv $bam_dir $outputdir
$commandsmax_se $tsv $bam_dir $outputdir
cd $outputdir
for ID in {1..13} {15..26}
do
    sample=$( cat $list | sed "${ID}q;d" )
    grep "$sample" filterreadcommands.txt  > filterreadcommands_${sample}.txt
done

grep "Intestine_rep2" filterreadcommands_se.txt | sed "s#/bar/hlee#home/hyungjoo.lee#g" > filterreadcommands_Intestine_rep2.txt

mv filterreadcommands_Brain_rep1.txt filterreadcommands_Brain_rep1.txt.bak
grep -v "embryonic" filterreadcommands_Brain_rep1.txt.bak > filterreadcommands_Brain_rep1.txt

mv filterreadcommands_Brain_rep2.txt filterreadcommands_Brain_rep2.txt.bak
grep -v "embryonic" filterreadcommands_Brain_rep2.txt.bak > filterreadcommands_Brain_rep2.txt


for f in filterreadcommands_*.txt;do sed -i "s#(1_of_many)#\\\(1_of_many\\\)#g" $f; done

### Calculation is done on the HTCF
module load samtools
for f in filterreadcommands_*.txt
do
	echo $f
	while read -r line
	do
		eval $line &
	done < "$f"
done

wait


### (5) find
find ./filterreadstats -maxdepth 1 -name "*.stats" -type f -print0 | 
    xargs -0 -n128 -P1 grep e | sed 's/.stats\:/.stats\t/g' > filter_read_stats.txt



