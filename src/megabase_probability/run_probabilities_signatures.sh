#!/usr/bin/env bash

bash src/megabase_probability/get_hg19_windows.sh
python src/megabase_probability/exonic_regions.py

python src/megabase_probability/generate_megabase_and_consequences_files.py

# calculate the probability at megabase level. This excludes CGC sites
python src/megabase_probability/probability_spectra_MB.py

# calculate probabilities in coding, CGC
python src/megabase_probability/probability_mutations_affecting_coding.py

