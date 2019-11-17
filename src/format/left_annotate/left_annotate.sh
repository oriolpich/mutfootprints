#!/usr/bin/env bash
set -e

# Left-align all mutations

input=$1
output_dir=$2

folder=$(basename $(dirname ${input}))
file=$(basename ${input})
name="${file/_post_processed.vcf.gz/}"

mkdir -p ${output_dir}/${folder}

bcftools norm --check-ref -x \
    --fasta-ref data/bcftools/Homo_sapiens.GRCh37.chroms.fa \
    ${input} \
    2> ${output_dir}/${folder}/${name}.la.stats.txt \
    | gzip 1> ${output_dir}/${folder}/${name}.la.vcf.gz
