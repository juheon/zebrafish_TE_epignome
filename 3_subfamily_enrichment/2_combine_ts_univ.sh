#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
##SBATCH --array=1-11
#SBATCH --job-name=lor_subfam

# SOFTWARE
module load bedtools/2.27.1


# GENOME FEATURES
subfam_bed=$(ls ../annotation/TEsubfamily/*bed.gz)
genome_size=../annotation/danRer10.chrom.sizes



# INPUT DATA
dir_in=`pwd`
data_prom=${dir_in}/subfam_enriched_active_promoter_[bchiklmstu]*.txt
data_weak=${dir_in}/subfam_enriched_weak_promoter_[bchiklmstu]*.txt
data_enhc=${dir_in}/subfam_enriched_enhancer_[bchiklmstu]*.txt

data_prox=${dir_in}/subfam_enriched_proximal_ATAC_[bchiklmstu]*.txt
data_dist=${dir_in}/subfam_enriched_distal_ATAC_[bchiklmstu]*.txt

data_umr=${dir_in}/subfam_enriched_UMR_[bchiklmstu]*.txt
data_lmr=${dir_in}/subfam_enriched_LMR_[bchiklmstu]*.txt


# OUTPUT
dir_out=`pwd`

out_prom=${dir_out}/subfam_enrich_ts_bp_filtered_active_promoter.txt
out_weak=${dir_out}/subfam_enrich_ts_bp_filtered_weak_promoter.txt
out_enhc=${dir_out}/subfam_enrich_ts_bp_filtered_active_enhancer.txt
out_hete=${dir_out}/subfam_enrich_ts_bp_filtered_heterochromatin.txt

out_prox=${dir_out}/subfam_enrich_ts_bp_filtered_proximal_ATAC.txt
out_dist=${dir_out}/subfam_enrich_ts_bp_filtered_distal_ATAC.txt

out_umr=${dir_out}/subfam_enrich_ts_bp_filtered_UMR.txt
out_lmr=${dir_out}/subfam_enrich_ts_bp_filtered_LMR.txt


# COMMANDS
### Filtering
### LOR > 2
cat $data_prom |
    grep -w -f <(cat $data_prom | awk '$1!="NA" && $1>2' | cut -f2 | sort -u ) - > $out_prom

cat $data_weak |
    grep -w -f <(cat $data_weak | awk '$1!="NA" && $1>2' | cut -f2 | sort -u ) - > $out_weak

cat $data_enhc |
    grep -w -f <(cat $data_enhc | awk '$1!="NA" && $1>2' | cut -f2 | sort -u ) - > $out_enhc

cat $data_hete |
    grep -w -f <(cat $data_hete | awk '$1!="NA" && $1>2' | cut -f2 | sort -u ) - > $out_hete

cat $data_prox |
    grep -w -f <(cat $data_prox | awk '$1!="NA" && $1>2' | cut -f2 | sort -u ) - > $out_prox

cat $data_dist |
    grep -w -f <(cat $data_dist | awk '$1!="NA" && $1>2' | cut -f2 | sort -u ) - > $out_dist

cat $data_umr |
    grep -w -f <(cat $data_umr | awk '$1!="NA" && $1>2' | cut -f2 | sort -u ) - > $out_umr

cat $data_lmr |
    grep -w -f <(cat $data_lmr | awk '$1!="NA" && $1>2' | cut -f2 | sort -u ) - > $out_lmr


cat $data_tad |
    grep -w -f <(cat $data_tad | awk '$1!="NA" && $1>2' | cut -f2 | sort -u ) - > $out_tad

cat $data_loop |
    grep -w -f <(cat $data_loop | awk '$1!="NA" && $1>2' | cut -f2 | sort -u ) - > $out_loop



