#!/bin/bash

python cotreatment.py run-create-random-matrix-treatments --deconstruction_folder ./data/deconstruction/breast --treatment_matrix_folder ./data/treatment_matrix/breast --symmetric
python cotreatment.py run-create-random-matrix-treatments --deconstruction_folder ./data/deconstruction/breast --treatment_matrix_folder ./data/treatment_matrix/breast
python cotreatment.py run-create-random-matrix-treatments --deconstruction_folder ./data/deconstruction/breastlungcolon --treatment_matrix_folder ./data/treatment_matrix/breastlungcolon --symmetric
python cotreatment.py run-create-random-matrix-treatments --deconstruction_folder ./data/deconstruction/breastlungcolon --treatment_matrix_folder ./data/treatment_matrix/breastlungcolon

