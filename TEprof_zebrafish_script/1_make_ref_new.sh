#!/bin/bash
# Author: Hyung Joo Lee


### SOFTWARE
ml htslib/1.3.1
gencode2dic=../TEprof_zebrafish/bin/genecode_to_dic.py
gencode_introns=../TEprof_zebrafish/bin/genecode_introns.py


### INPUT
in_gtf=../annotation/Danio_rerio.GRCz10.91.gtf
bed_rmsk=../annotation/repeats.bed.gz
in_rmsk=../annotation/rmsk.txt.gz


### OUTPUT
outdir=../output_teprof_zebrafish
out_gtf=$outdir/Danio_rerio.GRCz10.91.teprof_zebrafish.gtf
out_bed_rmsk=$outdir/rmsk.bed.gz
out_rmsk=$outdir/repeatmasker_danRer10_description_uniq.lst
out_intron_plus=$outdir/Danio_rerio.GRCz10.91.gtf_introns_plus.bed.gz
out_intron_minus=$outdir/Danio_rerio.GRCz10.91.gtf_introns_minus.bed.gz


### COMMANDS
# 1. Assemble transcripts and make GTF file
## transcript_id should be printed at the end of the row as a separate field.
## transcript_id is usually the thrid attribute in the Ensembl GTF file, but just in case...
## Some gene names have space in it "GENE NAME (1 of many)". Replace space character with _.
## Transcripts in the Unplaced scaffold are omiited for the downstream pipeline... (total 932 transcripts)
cat $in_gtf  |
    awk '$3=="transcript" || $3=="exon" || $3=="start_codon"' |
    sed "s/ (1 of /_(1_of_/g" | sed "/^chrUn/d" |
    awk -F "; " '{printf "%s",$0; for (i=2;i<=NF;i++) { if($i~/^transcript_id /) {printf "\t%s\n",$2}}}' > $out_gtf

# 2. Run custom python script
## Custom pyhon script modified
## COMMETED OUT: if (('appris_principal' in temp[-2])):
## Zebrafish Ensembl does not have appris_principal attribute.
## Changed the indentation accordingly.
$gencode2dic $out_gtf
sed -i "s/_number \"\([0-9]*\)\";/_number \1;/g" genecode_plus.dic
sed -i "s/_number \"\([0-9]*\)\";/_number \1;/g" genecode_minus.dic
mv genecode_plus.dic $odir
mv genecode_minus.dic $odir

# 3. Repeatmasker   ### 3,565,006 individual fragments / copies
# Downloaded from the UCSC genome browser
# TE in the Unplaced scafoold are excluded for the downstream pipeline... 
zcat $bed_rmsk | sed "/^chrUn/d" | sort -k1,1 -k2,2n | bgzip -c > $out_bed_rmsk     # 3,475,284 elements (loci)
tabix -p bed $out_bed_rmsk

# 4. Repeatmasker mapping
zcat $in_rmsk | awk '$6!~/^chrUn/' |  cut -f11-13 | sort -u > $out_rmsk   # 11,260 --> 11,137

# 5. Intron annotation
$gencode_introns $out_gtf
sort -k1,1 -k2,2n -k3,3n $out_gtf"_introns_plus"  | sed "s/_number \"\([0-9]*\)\";/_number \1;/g" | bgzip -c > $out_intron_plus
sort -k1,1 -k2,2n -k3,3n $out_gtf"_introns_minus" | sed "s/_number \"\([0-9]*\)\";/_number \1;/g" | bgzip -c > $out_intron_minus
tabix -p bed $out_intron_plus
tabix -p bed $out_intron_minus
rm $out_intron_plus $out_intron_minus


 

