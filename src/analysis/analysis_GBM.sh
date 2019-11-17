#!/usr/bin/env bash

mkdir -p data/GBM_Wang/

cd data/GBM_Wang/

wget https://www.nature.com/ng/journal/v48/n7/extref/ng.3590-S4.xlsx
wget https://media.nature.com/original/nature-assets/ng/journal/v48/n7/extref/ng.3590-S12.xlsx


cd ../../

python src/format/format_gbm_wang.py

Rscript src/signatures/deconstructsigs_sender.R  data/GBM_Wang/gbm_muts.to_deconstruct.tsv data/GBM_Wang/gbm_muts.to_deconstruct.out.tsv exome2genome
python src/analysis/GBM_process.py