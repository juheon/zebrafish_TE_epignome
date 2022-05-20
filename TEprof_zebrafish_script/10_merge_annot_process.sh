#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
#ml htslib/1.3.1
#ml python2/2.7

ml R/3.6.1

r_script=`pwd`/../TEprof_zebrafish/bin/mergeAnnotationProcess.R
argument=./arguments.txt
inputdir1=`pwd`/../output_teprof_zebrafish/08_annot_merged_gtf
inputdir2=`pwd`/../output_teprof_zebrafish/09_cal_trx_level_tpm
outputdir=`pwd`/../output_teprof_zebrafish/10_process_output
py_script=`pwd`/../TEprof_zebrafish/bin/stringtieExpressionFrac.py

log=$outputdir/mergeAnnotationProcess.R.out

istat=`pwd`/../output_teprof_zebrafish/06_filter_candi_read_info/Step6.RData

cd $inputdir1
$r_script -i $istat -o $outputdir &>$log

cd $inputdir2
## INTRON coverage info
find . -name "*i_data.ctab" > ctab_i.txt

cat ctab_i.txt | 
    while read ID ; do fileid=$(echo "$ID" | awk -F "/" '{print $2}'); cat <(printf 'chr\tstrand\tstart\tend\t'${fileid/_stats/}'\n') <(grep -F -f $outputdir/candidate_introns.txt $ID | awk -F'\t' '{ print $2"\t"$3"\t"$4"\t"$5"\t"$6 }') > ${ID}_cand ; done 

paste */*i_data.ctab_cand | awk '{ printf("%s\t%s\t%s\t%s", $1,$2,$3,$4); for (i=5;i<=NF;i+=5) {printf("\t%s",$i)} print "" }' > table_i_all

### transcript-level expression information

for f in ./*_stats/t_data.ctab
do
    $py_script $f &
done

wait

## Aggregate the stats
ls ./*stats/t_data.ctab_frac_tot > ctab_frac_tot_files.txt
ls ./*stats/t_data.ctab_tpm > ctab_tpm_files.txt

cat <(echo "TranscriptID") <(find . -name "*ctab_frac_tot" | head -1 | while read file ; do sort $file | awk '{print $1}' ; done;) > table_frac_tot
cat ctab_frac_tot_files.txt | while read file ; do fileid=$(echo "$file" | awk -F "/" '{print $2}') ; paste -d'\t' <(cat table_frac_tot) <(cat <(echo ${fileid/_stats/}) <(sort $file | awk '{print $2}')) > table_frac_tot_temp; mv -f table_frac_tot_temp table_frac_tot; done

cat <(echo "TranscriptID") <(find . -name "*ctab_tpm" | head -1 | while read file ; do sort $file | awk '{print $1}' ; done;) > table_tpm
cat ctab_tpm_files.txt | while read file ; do fileid=$(echo "$file" | awk -F "/" '{print $2}') ; paste -d'\t' <(cat table_tpm) <(cat <(echo ${fileid/_stats/}) <(sort $file | awk '{print $2}')) > table_tpm_temp; mv -f table_tpm_temp table_tpm; done

cat <(head -1 table_frac_tot) <(grep -Ff $outputdir/candidate_names.txt table_frac_tot) > table_frac_tot_cand
cat <(head -1 table_tpm) <(grep -Ff $outputdir/candidate_names.txt table_tpm) > table_tpm_cand

mv table_* $outputdir
mv ctab_* $outputdir

