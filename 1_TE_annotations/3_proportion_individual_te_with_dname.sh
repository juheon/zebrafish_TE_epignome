#!/bin/bash
# Author: Hyung Joo Lee

#SBATCH --mem=2G
#SBATCH --array=1-11
##SBATCH --workdir=/scratch/twlab/hlee/zf_te
#SBATCH --job-name=prop_te_dname

# SOFTWARE
module load bedtools/2.27.1
ID=$SLURM_ARRAY_TASK_ID


# TE classes
dir_te=/scratch/twlab/hlee/genomes/danRer10/rmsk
bed_dna=${dir_te}/danRer10.DNA.bed.gz
bed_ltr=${dir_te}/danRer10.LTR.bed.gz
bed_line=${dir_te}/danRer10.LINE.bed.gz
bed_sine=${dir_te}/danRer10.SINE.bed.gz
bed_rc=${dir_te}/danRer10.RC.bed.gz
bed_te=${dir_te}/danRer10.TE.bed.gz


# INPUT DATA
list=tissues-e.txt
tissue=$( cat $list | sed "${ID}q;d" )


te_me=2_dynamic/danRer10.TE_wCG.DNAme_${tissue}_cov5.bed.gz


# OUTPUT
dir_out=2_dynamic
out1=${dir_out}/fraction_dname_states_${tissue}.txt
out2=${dir_out}/fraction_dname_levels_${tissue}.txt


# COMMANDS
#cnt_te=$( zcat $bed_te | wc -l )

#umr_dna=$(intersectBed -sorted -u -a $bed_dna -b $umr | wc -l )
#umr_ltr=$(intersectBed -sorted -u -a $bed_ltr -b $umr | wc -l )
#umr_line=$(intersectBed -sorted -u -a $bed_line -b $umr | wc -l )
#umr_sine=$(intersectBed -sorted -u -a $bed_sine -b $umr | wc -l )
#umr_te=$(( umr_dna + umr_ltr + umr_line + umr_sine ))

#awk -v"a=$umr_te" -v"b=$cnt_te" -v"n=$tissue" 'BEGIN{printf("%.3g\t%d\t%d\tUMR\t%s\n", a/b,a,b,n)}' > $out1

#lmr_dna=$(intersectBed -sorted -u -a $bed_dna -b $lmr | wc -l )
#lmr_ltr=$(intersectBed -sorted -u -a $bed_ltr -b $lmr | wc -l )
#lmr_line=$(intersectBed -sorted -u -a $bed_line -b $lmr | wc -l )
#lmr_sine=$(intersectBed -sorted -u -a $bed_sine -b $lmr | wc -l )
#lmr_te=$(( lmr_dna + lmr_ltr + lmr_line + lmr_sine ))

#awk -v"a=$lmr_te" -v"b=$cnt_te" -v"n=$tissue" 'BEGIN{printf("%.3g\t%d\t%d\tLMR\t%s\n", a/b,a,b,n)}' >> $out1

zcat $te_me |
    awk -v n=$tissue 'BEGIN{na=0; lo=0; me=0; hi=0}
            $5=="NA" { na++ } $5!="NA" && $5<0.25 { lo++ }
            $5!="NA" && $5>=0.25 && $5<0.75 { me++ } $5!="NA" && $5>=0.75 { hi++ }
         END{ printf("%3g\t%d\t%d\tMissing\t%s\n", na/NR,na,NR,n );
              printf("%3g\t%d\t%d\tLow\t%s\n", lo/NR,lo,NR,n );
              printf("%3g\t%d\t%d\tMedium\t%s\n", me/NR,me,NR,n );
              printf("%3g\t%d\t%d\tHigh\t%s\n", hi/NR,hi,NR,n ) }' > $out2


