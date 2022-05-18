#!/bin/bash
# Author: Hyung Joo Lee

# SOFTWARE
module load bedtools/2.27.1

# INPUT
rmsk=rmsk.txt.gz

cpg=CpG_sites.bed.gz

# OUTPUT
dna=danRer10.DNA.bed.gz		## DNA + DNA?
ltr=danRer10.LTR.bed.gz
line=danRer10.LINE.bed.gz
sine=danRer10.SINE.bed.gz
rc=danRer10.RC.bed.gz		## RC
te=danRer10.TE.bed.gz		## all 

te_cg=danRer10.TE_wCG.bed.gz
te_frag=danRer10.TE_frag.bed.gz

out=danRer10.CpGs_in_TE.txt

# COMMANDS
# REPEAT by CLASS
zcat $rmsk | 
	awk -F"\t" -vOFS="\t" '$12~/^DNA/ {print $6,$7,$8}' | 
	sort -k1,1 -k2,2n |
	mergeBed -i stdin | gzip -nc > $dna

zcat $rmsk | 
	awk -F"\t" -vOFS="\t" '$12=="LTR" {print $6,$7,$8}' | 
	sort -k1,1 -k2,2n |
	mergeBed -i stdin | gzip -nc > $ltr

zcat $rmsk | 
	awk -F"\t" -vOFS="\t" '$12=="LINE" {print $6,$7,$8}' | 
	sort -k1,1 -k2,2n |
	mergeBed -i stdin | gzip -nc > $line

zcat $rmsk | 
	awk -F"\t" -vOFS="\t" '$12=="SINE" {print $6,$7,$8}' | 
	sort -k1,1 -k2,2n |
	mergeBed -i stdin | gzip -nc > $sine

zcat $rmsk | 
	awk -F"\t" -vOFS="\t" '$12=="RC" {print $6,$7,$8}' | 
	sort -k1,1 -k2,2n |
	mergeBed -i stdin | gzip -nc > $rc

zcat $dna $ltr $line $sine $rc |
	sort -k1,1 -k2,2n |
	mergeBed -i stdin | gzip -nc > $te


### TE with CGs 1,941,161 fragments
zcat $rmsk |
	awk -F"\t" -vOFS="\t" '$12~/^DNA/ || $12=="LTR" || $12~/[LS]INE/ || $12=="RC" {print $6,$7,$8,$11,0,$10,$12}' |
	sed "s/?$//g" | sort -k1,1 -k2,2n |
	intersectBed -u -sorted -F 1 -a stdin -b $cpg | gzip -nc > $te_cg

### TE fragments
zcat $rmsk |
	awk -F"\t" -vOFS="\t" '$12~/^DNA/ || $12=="LTR" || $12~/[LS]INE/ || $12=="RC" {print $6,$7,$8,$11,0,$10,$12}' |
	sed "s/?$//g" | sort -k1,1 -k2,2n | gzip -nc > $te_frag


for bed in ${dir}/danRer10.{TE,DNA,LTR,LINE,SINE,RC}.bed.gz
do
    base=${bed##*/}
    base=${base#danRer10.}
    base=${base%.bed.gz}
    cnt=$( intersectBed -sorted -u -f 1 -a $cpg -b $bed | wc -l )
    echo -e "$base\t$cnt" >> $out
done

