#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=2,8
#SBATCH --job-name=hic


# SOFTWARE
module load bedtools/2.27.1

ID=$SLURM_ARRAY_TASK_ID
workdir=`pwd`


### INPUT
tissue_list=${workdir}/tissues-e.txt
dir_in=${workdir}/0_cres

motif=ctcf.sites.danRer10.bed

tissue=$( cat $tissue_list | sed "${ID}q;d" )

boundaries=${dir_in}/TADs/TAD_boundaries_${tissue}.bed.gz
anchors=${dir_in}/loop_anchors/loop_anchors_${tissue}.bed.gz

dist=${dir_in}/distalatacseq/${tissue}.distal.open.bed
prox=${dir_in}/proximalatacseq/${tissue}.proximal.open.bed

### OUTPUT
atac_ctcf=${dir_in}/${tissue}.ctcf_atac.bed.gz

#tad_out=${dir_in}/TADs/TAD_boundaries_${tissue}.ctcf_atac.bed.gz
#loop_out=${dir_in}/loop_anchors/loop_anchors_${tissue}.ctcf_atac.bed.gz
tad_out=${dir_in}/TADs/TAD_boundaries_${tissue}.atac_ctcf.bed.gz
loop_out=${dir_in}/loop_anchors/loop_anchors_${tissue}.atac_ctcf.bed.gz


### COMMANDS
# TAD
cat $dist $prox | sort -k1,1 -k2,2n |
	intersectBed -sorted -wo -F 1 -a stdin -b $motif |
	groupBy -c 9 -o max | 
	intersectBed -sorted -wo -F 1 -a stdin -b $motif |
	awk -F"\t" -vOFS="\t" '$4==$9 {print $1,$2,$3,$5,$6,$7,$9}' |  gzip -nc > $atac_ctcf

#intersectBed -sorted -wo -F 0.2 -a $boundaries -b $atac_ctcf |
#	groupBy -c 10 -o max |
#	intersectBed -sorted -wo -a stdin -b $atac_ctcf |
#	awk -F"\t" -vOFS="\t" '$4==$11 {print $5,$6,$7,$8,$9,$10,$11}' | gzip -nc > $tad_out

intersectBed -sorted -u -a $atac_ctcf -b $boundaries | cut -f1-3 | gzip -nc > $tad_out

# loops
#intersectBed -sorted -wo -F 0.2 -a $anchors -b $atac_ctcf |
#	groupBy -c 10 -o max |
#	intersectBed -sorted -wo -a stdin -b $atac_ctcf |
#	awk -F"\t" -vOFS="\t" '$4==$11 {print $5,$6,$7,$8,$9,$10,$11}' | gzip -nc > $loop_out

intersectBed -sorted -u -a $atac_ctcf -b $anchors | cut -f1-3 | gzip -nc > $loop_out
