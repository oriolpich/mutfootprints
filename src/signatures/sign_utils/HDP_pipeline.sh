#!/usr/bin/env bash

# create input
python src/signatures/sign_utils/get_HDP_input.python
mkdir data/hartwig/signatures/HDP
# here we are showing how it would run, we strongly recommend running this in parallel
for i in {1..14}
do
  Rscript src/signatures/sign_utils/run_HDP.R $i
done

# merge results
Rscript src/signatures/sign_utils/results_HDP.R

# get the similarity between new signatures and SignatureAnalyzer
python src/signatures/sign_utils/similarity_HDP.py