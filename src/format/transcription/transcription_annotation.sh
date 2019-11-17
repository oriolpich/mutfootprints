#!/usr/bin/env bash

# Path of the data
PATH_FILE="data/asymmetry_files/"
mkdir -p ${PATH_FILE}
cd $PATH_FILE

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.annotation.gff3.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.annotation.gtf.gz
zgrep protein_coding gencode.v29lift37.annotation.gff3.gz | grep -w gene | cut -f1,4,5,7 | sort -k1,1 -k2,2n | awk '{OFS="\t"}{print $1, $2, $3, "1", "1", $4}' | grep "+" | bedtools merge -i stdin -s -c 6 -o distinct | gzip > positive_strand_genes.bed.gz
zgrep protein_coding gencode.v29lift37.annotation.gff3.gz | grep -w gene | cut -f1,4,5,7 | sort -k1,1 -k2,2n | awk '{OFS="\t"}{print $1, $2, $3, "1", "1", $4}' | grep "-" | bedtools merge -i stdin -s -c 6 -o distinct | gzip > negative_strand_genes.bed.gz

intersectBed -a positive_strand_genes.bed.gz -b negative_strand_genes.bed.gz -v | gzip > positive_strand_nooverlapp.bed.gz
intersectBed -a negative_strand_genes.bed.gz -b positive_strand_genes.bed.gz -v | gzip > negative_strand_nooverlapp.bed.gz

zcat negative_strand_nooverlapp.bed.gz positive_strand_nooverlapp.bed.gz | sort -k1,1 -k2,2n | gzip > strand_coordinates.bed.gz
