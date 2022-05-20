grep "srl" candidate_transcripts.gff3 | fastaFromBed -split -s -fi /bar/genomes/danRer10/danRer10.fa -bed stdin > srl.fasta
