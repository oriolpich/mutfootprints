
# format mutations
python src/format/format_stjude.py

# run deconstructsigs
Rscript src/signatures/deconstructsigs_sender.R data/STJUDE/format/STJUDE.to_deconstruct.tsv data/STJUDE/format/STJUDE.to_deconstruct.out.tsv genome

# plot
python src/analysis/fitting_stjude.py
