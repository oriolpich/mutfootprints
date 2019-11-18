#!/bin/bash

echo "SBS9 breast"
python logistic_regression.py treatment-shuffle --sig SBS9 --deconstruction_folder ./data/deconstruction/breast

echo "SBS9 breastlungcolon"
python logistic_regression.py treatment-shuffle --sig SBS9 --deconstruction_folder ./data/deconstruction/breastlungcolon

echo "SBS31 breast"
python logistic_regression.py treatment-shuffle --sig SBS31 --deconstruction_folder ./data/deconstruction/breast

echo "SBS31 breastlungcolon"
python logistic_regression.py treatment-shuffle --sig SBS31 --deconstruction_folder ./data/deconstruction/breastlungcolon
