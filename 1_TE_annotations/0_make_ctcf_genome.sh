#!/bin/bash
# Author: Hyung Joo Lee

# SOFTWARE
module load homer/4.8

motif=ctcf_zf.motif

scanMotifGenomeWide.pl $motif danRer10 -bed -keepAll > ctcf_zebrafish.sites.danRer10.bed

ctcf_zf_motif=CTCF_zebrafish.sites.bed.gz
cat ctcf_zebrafish.sites.danRer10.bed | awk -F"\t" -vOFS="\t" '{$2--; print }' |  sort -k1,1 -k2,2n | gzip -nc > $ctcf_zf_motif


zcat GSE133437_CTCF_24hpf_combined_reps_all_peaks.bed.gz | sed 1d | cut -f1-3 | sort -k1,1 -k2,2n | gzip -nc > ChIP_CTCF_24hpf_peaks.bed.gz
zcat GSE133437_CTCF_24hpf_combined_reps_all_peaks.bed.gz | sed 1d | awk -F"\t" -vOFS="\t" '{print $1,$2+$10,$2+$10+1}' | sort -k1,1 -k2,2n | gzip -nc > ChIP_CTCF_24hpf_summits.bed.gz

