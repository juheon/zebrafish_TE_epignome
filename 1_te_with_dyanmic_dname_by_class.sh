#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=annot_feat_epi

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID


# TE classes
dir_te=/scratch/twlab/hlee/genomes/danRer10/rmsk
bed_te_cg=${dir_te}/danRer10.TE_wCG.bed.gz

# INPUT DATA
dir_in=0_cres
te_me=1_annot/danRer10.TE_wCG.DNAme_*_cov5.bed.gz


# OUTPUT
dir_out=2_dynamic
out=${dir_out}/te_dynamic_dname_by_class.txt


# COMMANDS
intersectBed -sorted -wo -f 1 -F 1 -a $bed_te_cg -b $te_me |
    groupBy -g 1,2,3,7 -c 13 -o collapse | sed "s/,/\t/g" |
    awk -F"\t" ' { na=lo=me=hi=0;
        for(i=5; i<14; i++) {
            if( $i=="NA" ) { na=1 }
            if( $i!="NA" && $i<0.25 ) { lo=1 }
            if( $i!="NA" && $i>=0.25 && $i<0.75 ) { me=1 }
            if( $i!="NA" && $i>=0.75 ) { hi=1 } } }
        $4=="DNA" { dna[na+lo+me+hi]++ }
        $4=="LTR" { ltr[na+lo+me+hi]++ }
        $4=="LINE" { line[na+lo+me+hi]++ }
        $4=="SINE" { sine[na+lo+me+hi]++ }
    END{ OFS="\t"; print 1,2,3,4;
        print "DNA", dna[1], dna[2], dna[3], dna[4];
        print "LTR", ltr[1], ltr[2], ltr[3], ltr[4];
        print "LINE", line[1], line[2], line[3], line[4];
        print "SINE", sine[1], sine[2], sine[3], sine[4] } ' > $out

