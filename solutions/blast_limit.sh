#!/bin/bash

# Prep subject sequence
wget -N --no-check-certificate "https://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_sm.toplevel.fa.gz"

gunzip -kf Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_sm.toplevel.fa.gz

makeblastdb -dbtype nucl -in Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_sm.toplevel.fa -input_type fasta -parse_seqids -title ecoli

# Prep query sequence
wget -N --no-check-certificate "https://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz"

gunzip -kf Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz

perl -pe '/^>/ ? print "\n" : chomp' \
Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa \
> Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.unwrap.fa

sed 1d Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.unwrap.fa \
| head -6 | tail -2 > mysequence.fa

# blast with unmodified sequence
blastn -evalue 0.001 \
-db Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_sm.toplevel.fa \
-query mysequence.fa \
-outfmt 6

# initial run
for i in $(seq 10 10 800) ; do

  echo No point mutations: $i

  msbar -sequence mysequence.fa -count $i -point 4 -block 0 -codon 0 \
  -outseq /dev/stdout 2>/dev/null > mysequence_mutated.fa

  blastn -evalue 0.001 \
  -db Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_sm.toplevel.fa \
  -query mysequence_mutated.fa \
  -outfmt 6
done




# now want to run replications at each value of i
for i in $(seq 0 50 700) ; do

  RES=$(for j in {1..100} ; do

    msbar -sequence mysequence.fa -count $i -point 4 -block 0 -codon 0 \
    -outseq /dev/stdout 2>/dev/null > mysequence_mutated.fa

    blastn -evalue 0.001 \
    -db Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_sm.toplevel.fa \
    -query mysequence_mutated.fa \
    -outfmt 6 | grep -w AAC73114 | head -1 | wc -l

  done | numaverage)

  echo $i $RES | tr ' ' '\t'

done | tee blastres.tsv



# follow up in area of interest
for i in $(seq 200 10 300) ; do

  RES=$(for j in {1..100} ; do

    msbar -sequence mysequence.fa -count $i -point 4 -block 0 -codon 0 \
    -outseq /dev/stdout 2>/dev/null > mysequence_mutated.fa

    blastn -evalue 0.001 \
    -db Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_sm.toplevel.fa \
    -query mysequence_mutated.fa \
    -outfmt 6 | grep -w AAC73114 | head -1 | wc -l

  done | numaverage)

  echo $i $RES | tr ' ' '\t'

done | tee blastres2.tsv


