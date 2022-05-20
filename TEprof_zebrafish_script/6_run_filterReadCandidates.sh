#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
#ml htslib/1.3.1
#ml python2/2.7

ml R/3.6.1

r_script=`pwd`/../TEprof_zebrafish/bin/filterReadCandidates.R

inputdir1=`pwd`/../output_teprof_zebrafish/04_aggregate_annotation_process_output
inputdir2=`pwd`/../output_teprof_zebrafish/05_cal_read_info
outputdir=`pwd`/../output_teprof_zebrafish/06_filter_candi_read_info
mkdir $outputdir

read_info=$inputdir2/filter_read_stats.txt

argument=./arguments.txt



### INPUT # 26 annotated all_c
#c_files=$( ls ${workdir}/*_rep1.gtf_annotated_filtered_test_all_c )

log=$outputdir/filterReadCandidates.R.out
cd $inputdir2
step4data=$inputdir1/Step4.RData
$r_script -i $step4data -t $read_info -o $outputdir &>$log
wait



## Check the transcripts by sample from the step 5
## output from the above: read_filtered_candidates.tsv
## 

candidates=$outputdir/read_filtered_candidates.tsv
samples=$inputdir1/samples.txt

## output
cand_read_info=$outputdir/read_filtered_candidates_read_info_by_sample.txt


numCand=$( sed 1d $candidates | wc -l )
## 845 --> 946




for num in {1..946}
do
    cand=$( sed 1d $candidates | sed "${num}q;d" | cut -f1 )
    cnt=$( sed 1d $candidates | sed "${num}q;d" | cut -f18 )
    grep -e "$cand" $read_info | sed "s#^./filterreadstats/##g; s/--/\t/g; s/.stats\t/\t/g" |
        sort -k3,3nr -k4,4nr | head -${cnt}
done > $cand_read_info


cut $f2 $cand_read_info | sort | uniq -c



