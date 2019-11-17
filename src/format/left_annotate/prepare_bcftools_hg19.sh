#!/usr/bin/env bash

# path where the final file will be stored. Substitute it with you own path
PATH_FILE="data/bcftools"
mkdir -p ${PATH_FILE}

# download file from ensembl
wget -O ${PATH_FILE}/Homo_sapiens.GRCh37.dna_sm.toplevel.fa.gz ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_sm.toplevel.fa.gz

# extract it
zcat ${PATH_FILE}/Homo_sapiens.GRCh37.dna_sm.toplevel.fa.gz > ${PATH_FILE}/Homo_sapiens.GRCh37.dna_sm.toplevel.fa

# create index
samtools faidx ${PATH_FILE}/Homo_sapiens.GRCh37.dna_sm.toplevel.fa

# merge chromosomes
samtools faidx ${PATH_FILE}/Homo_sapiens.GRCh37.dna_sm.toplevel.fa 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT> ${PATH_FILE}/Homo_sapiens.GRCh37.chroms.fa

# create index for bcftools
samtools faidx ${PATH_FILE}/Homo_sapiens.GRCh37.chroms.fa
