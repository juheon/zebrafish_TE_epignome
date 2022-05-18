#!/bin/bash
# Author: Yiran Hou

#### Age comparison
# Intersect active regions with age estimation
for i in danRer10.TEs_11_tissues.*.bed; do OUT="danRer10.RM.JCdist.${i#*TEs_}";echo $OUT; bedtools intersect -c -a danRer10.TE.copy.age.txt -b $i > $OUT;done 
for i in *bed;do OUT="${i#*tissues.}.ct";cut -f10 $i > $OUT;done
awk '$1==0{print "Others"}$1!=0{print "ActiveEnhancer"}' active_enh.bed.ct > active_enh.bed.ct.nm 
awk '$1==0{print "Others"}$1!=0{print "ActivePromoter"}' active_pro.bed.ct > active_pro.bed.ct.nm
awk '$1==0{print "Others"}$1!=0{print "ATAC"}' ATAC.bed.ct > ATAC.bed.ct.nm
awk '$1==0{print "Others"}$1!=0{print "CRE"}' CRE.bed.ct > CRE.bed.ct.nm
awk '$1==0{print "Others"}$1!=0{print "DistalATAC"}' distal_ATAC.bed.ct > distal_ATAC.bed.ct.nm
awk '$1==0{print "Others"}$1!=0{print "Heterochromatin"}' heterochromatin.bed.ct > heterochromatin.bed.ct.nm
awk '$1==0{print "Others"}$1!=0{print "LMR"}' LMR.bed.ct > LMR.bed.ct.nm
awk '$1==0{print "Others"}$1!=0{print "ProximalATAC"}' proximal_ATAC.bed.ct > proximal_ATAC.bed.ct.nm
awk '$1==0{print "Others"}$1!=0{print "UMR"}' UMR.bed.ct > UMR.bed.ct.nm
awk '$1==0{print "Others"}$1!=0{print "UMR_LMR"}' UMR_LMR.bed.ct > UMR_LMR.bed.ct.nm
awk '$1==0{print "Others"}$1!=0{print "WeakPromoter"}' weak_pro.bed.ct > weak_pro.bed.ct.nm
## Paste all annotations together
paste danRer10.RM.JCdist active_enh.bed.ct.nm CRE.bed.ct.nm LMR.bed.ct.nm UMR_LMR.bed.ct.nm active_pro.bed.ct.nm distal_ATAC.bed.ct.nm proximal_ATAC.bed.ct.nm weak_pro.bed.ct.nm ATAC.bed.ct.nm heterochromatin.bed.ct.nm UMR.bed.ct.nm > danRer10.RM.JCdist.11_tissues.allAnn
## Split TE family and subfamily 
awk '{OFS="\t";split($4,a,"#");print $1,$2,$3,$4,a[1],a[2],$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' danRer10.RM.JCdist.11_tissues.allAnn > danRer10.RM.JCdist.11_tissues.AllAnnd
## Clean up to include only the major five classes
awk '{OFS="\t";split($6,a,"/"); print $1,$2,$3,$4,$5,$6,a[1],$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' danRer10.RM.JCdist.11_tissues.AllAnnd | awk '($7=="DNA" || $7=="LTR" || $7=="LINE" || $7=="SINE" || $7=="RC"){OFS="\t"; print $0}' | sort -k1,1 -k2,2n | uniq > danRer10.RM.JCdist.11_tissues.cleaned.AllAnnd
## Add testis TETSS as a column in metadata with age 
cut -f1-3 TETSS_wATAC.TestisUniq.bed | uniq > ATAC_filtered_Testis_TSS.bed
bedtools intersect -c -a danRer10.RM.JCdist.11_tissues.cleaned.AllAnnd -b ATAC_filtered_Testis_TSS.bed > danRer10.RM.JCdist.11_tissues.filtered.TestisTSS | awk '{OFS="\t";for(i=1;i<=23;i++) printf("%s%s", $i, OFS)}$24==0{print "Others"}$24!=0{print "Testis_TETSS"}' | awk '$11<=0.5{print $0}' | sed 's/-0.0/0.0/g' > danRer10.RM.JCdist.11_tissues.filtered.TestisTSS.sub05.ann

## Extract all active states
awk '$13=="ActiveEnhancer"{OFS="\t"; print $7,$12,$13}' danRer10.RM.JCdist.11_tissues.cleaned.TestisTSS.sub05.ann > nm.act_enh &
awk '$17=="ActivePromoter"{OFS="\t"; print $7,$12,$17}' danRer10.RM.JCdist.11_tissues.cleaned.TestisTSS.sub05.ann > nm.active_pro &
awk '$20=="WeakPromoter"{OFS="\t"; print $7,$12,$20}' danRer10.RM.JCdist.11_tissues.cleaned.TestisTSS.sub05.ann > nm.weak_pro &
awk '$19=="ProximalATAC"{OFS="\t"; print $7,$12,$19}' danRer10.RM.JCdist.11_tissues.cleaned.TestisTSS.sub05.ann > nm.pro_atac &
awk '$18=="DistalATAC"{OFS="\t"; print $7,$12,$18}' danRer10.RM.JCdist.11_tissues.cleaned.TestisTSS.sub05.ann > nm.dist_atac &
awk '$23=="UMR"{OFS="\t"; print $7,$12,$23}' danRer10.RM.JCdist.11_tissues.cleaned.TestisTSS.sub05.ann > nm.umr &
awk '$15=="LMR"{OFS="\t"; print $7,$12,$15}' danRer10.RM.JCdist.11_tissues.cleaned.TestisTSS.sub05.ann > nm.lmr & 
awk '($14=="Others" && $16=="Others" && $21=="Others"){OFS="\t"; print $7,$12,"Others"}' danRer10.RM.JCdist.11_tissues.cleaned.TestisTSS.sub05.ann > nm.inactive &
## Concatenate all annotations together
cat nm.act_enh nm.active_pro nm.weak_pro nm.pro_atac nm.dist_atac nm.umr nm.lmr nm.inactive > nm.anns

