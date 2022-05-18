#!/bin/bash
# Author: Hyung Joo Lee

# SOFTWARE
module load kentsrc/20170117
module load bedtools/2.27.1

# INPUT
ensembl=Danio_rerio.GRCz10.91.gtf.gz
chromSize=danRer10.chrom.sizes 

cgi_ucsc=cpgIslandExt.txt.gz

# TEMPORARY OUTPUT
genpred=features/GRCz10.91.annotation.genPred
tmp=features/GRCz10.91.annotation.unsrt.bed

# OUTPUT
bed12=features/GRCz10.91.annotation.bed.gz
exon=features/GRCz10.91.exon.bed.gz
intron=features/GRCz10.91.intron.bed.gz
utr5=features/GRCz10.91.5utr.bed.gz
utr3=features/GRCz10.91.3utr.bed.gz
cds=features/GRCz10.91.cds.bed.gz
tss2k=features/GRCz10.91.tss2k.bed.gz
cgi=features/danRer10.cgi.bed.gz
intergenic=features/GRCz10.91.intergenic.bed.gz

# COMMANDS
## initial prep step
gtfToGenePred $ensembl $genpred
genePredToBed $genpred $tmp
sort -k1,1 -k2,2n $tmp | gzip -nc > $bed12
rm $genpred $tmp


# CODE from http://onetipperday.sterding.com/2012/11/get-intron-utr-cds-from-bed12-format.html
# INTRON
zcat $bed12 |
	awk '{OFS="\t"; split($11,a,","); split($12,b,","); A=""; B=""; 
		for(i=1;i<length(a)-1;i++) {A=A""(b[i+1]-b[i]-a[i])","; B=B""(b[i]+a[i]-(b[1]+a[1]))",";} 
		if($10>1) print $1,$2+a[1], $3-a[length(a)-1], $4,$5,$6,$2+a[1], $3-a[length(a)-1],$9,$10-1,A,B;}' |
	bedtools bed12tobed6 -i stdin | sort -k1,1 -k2,2n | gzip -nc >$intron

# 5'UTR
zcat $bed12 |
	awk '{OFS="\t";split($11,blockSizes,","); split($12,blockStarts,","); blockCount=$10; A=""; B="";
		if($7==$8) next; N=0; 
		if($6=="+" && $2<$7) {start=$2; end=$7; for(i=1;i<=blockCount;i++) if(($2+blockStarts[i]+blockSizes[i])<=$7) {A=A""blockSizes[i]",";B=B""blockStarts[i]","; end=($2+blockStarts[i]+blockSizes[i]); N++;} else { if(($2+blockStarts[i])<$7) {A=A""($7-$2-blockStarts[i])",";B=B""blockStarts[i]","; N++; end=$7;} break; } print $1,start,end,$4,$5,$6,start,end,$9,N,A,B;} if($6=="-" && $8<$3) {start=$8;end=$3; for(i=1;i<=blockCount;i++) if(($2+blockStarts[i])>=$8) {if(start==0) {A=blockSizes[i];B=0; start=$2+blockStarts[i];} else {A=A","blockSizes[i];B=B","($2+blockStarts[i]-start);} N++;} else { if(($2+blockStarts[i]+blockSizes[i])>$8) { A=($2+blockStarts[i]+blockSizes[i]-$8);B=0; N++; start=$8;} if(($2+blockStarts[i]+blockSizes[i])==$8) start=0;} print $1,start,end,$4,$5,$6,start,end,$9,N,A,B;}}' |
	bedtools bed12tobed6 -i stdin | sort -k1,1 -k2,2n | gzip -nc >$utr5

# 3'UTR
zcat $bed12 |
	awk '{OFS="\t";split($11,blockSizes,","); split($12,blockStarts,","); blockCount=$10;A=""; B=""; if($7==$8) next;N=0;if($6=="-" && $2<$7) {start=$2;end=$7; for(i=1;i<=blockCount;i++) if(($2+blockStarts[i]+blockSizes[i])<=$7) {A=A""blockSizes[i]",";B=B""blockStarts[i]","; end=($2+blockStarts[i]+blockSizes[i]); N++;} else { if(($2+blockStarts[i])<$7) {A=A""($7-$2-blockStarts[i])",";B=B""blockStarts[i]","; N++; end=$7;} break; } print $1,start,end,$4,$5,$6,start,end,$9,N,A,B;} if($6=="+" && $8<$3) {start=$8;end=$3; for(i=1;i<=blockCount;i++) if(($2+blockStarts[i])>$8) { if(start==0) {A=blockSizes[i];B=0; start=$2+blockStarts[i];} else {A=A","blockSizes[i];B=B","($2+blockStarts[i]-start);} N++; } else { if(($2+blockStarts[i]+blockSizes[i])>$8) { A=($2+blockStarts[i]+blockSizes[i]-$8);B=0; N++; start=$8;} if(($2+blockStarts[i]+blockSizes[i])==$8) start=0;} print $1,start,end,$4,$5,$6,start,end,$9,N,A,B;}}' |
	bedtools bed12tobed6 -i stdin | sort -k1,1 -k2,2n | gzip -nc >$utr3

# EXON
zcat $ensembl |
	awk 'BEGIN{OFS="\t"} $3=="exon" { print $1,$4-1,$5,$14,"0",$7 }' |
	sed 's/\";//g; s/\"//g' |
	sort -k1,1 -k2,2n | uniq | gzip -nc > $exon

# CDS
zcat $ensembl |
	awk 'BEGIN{OFS="\t"} $3=="CDS" { print $1,$4-1,$5,$14,"0",$7 }' |
	sed 's/\";//g; s/\"//g' |
	sort -k1,1 -k2,2n | uniq | gzip -nc > $cds

# PROMOTERS (all genes ? or all transcript?)
zcat $ensembl |
	awk 'BEGIN{OFS="\t"} $3=="transcript" { if($7=="+") {$5=$4+500; $4=$4-2000} 
                                                if($7=="-") {$4=$5-500; $5=$5+2000} 
                                                if($4<0) {$4=0} print $1,$4,$5,$14,"0",$7 }' |
    sed 's/\";//g; s/\"//g' | sort -k1,1 -k2,2n |
    intersectBed -sorted -a stdin -b <( cat $chromSize | awk -vOFS="\t" '{print $1,0,$2}' | sort -k1,1 ) |
    uniq | gzip -nc > $promoter

# TSS +-1kb (all transcript)
zcat $ensembl |
	awk 'BEGIN{OFS="\t"} $3=="transcript" { if($7=="+") {$5=$4+1000; $4=$4-1000} 
                                                if($7=="-") {$4=$5-1000; $5=$5+1000} 
                                                if($4<0) {$4=0} print $1,$4,$5,$14,"0",$7 }' |
    sed 's/\";//g; s/\"//g' | sort -k1,1 -k2,2n |
    intersectBed -sorted -a stdin -b <( cat $chromSize | awk -vOFS="\t" '{print $1,0,$2}' | sort -k1,1 ) |
    uniq | sort -k1,1 -k2,2n | gzip -nc > $tss2k


# CGI
zcat $cgi_ucsc |
	awk -F"\t" 'BEGIN{OFS="\t"} { print $2,$3,$4,$5,"0","+"}' |
	sed 's/CpG: /CGI_/g' |
	sort -k1,1 -k2,2n | uniq | gzip -nc > $cgi

# INTERGENIC
chromSize=${dir}/danRer10.chrom.sizes
zcat $promoter $exon $intron $utr5 $utr3 $cds | sort -k1,1 -k2,2n |
    mergeBed -i stdin |
    subtractBed -sorted -a <( cat $chromSize | awk -vOFS="\t" '{print $1,0,$2}' | sort -k1,1 ) -b stdin |
    gzip -nc > $intergenic

