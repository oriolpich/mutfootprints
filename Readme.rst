
README
======

|

This folder contains the code used in `The Mutational Footprints of Cancer Therapies <LINK TO THE PAPER>`_.
If you use this code in a publication, please cite:

.. admonition:: Citation
   :class: note

   Oriol Pich, Ferran Muiños, Martijn Paul Lolkema, Neeltje Steeghs, Abel Gonzalez-Perez, Nuria Lopez-Bigas, The mutational footprints of cancer therapies, Nature Genetics, 2019 doi: `XX <LINK>`_

As part of this work we didn’t generate any original data. We re-used publicly available data described in specific sections of the methods.
The metastatic tumor cohort data (DR-024 version 2) is available from the `Hartwig Medical Foundation <(https://www.hartwigmedicalfoundation.nl/en>`_ for academic research upon request. Without the Hartwig Metastatic tumor cohort the results won't be reproducible.

|

|

Running this software
*********************

These analysis have been perform using software in Python, R, Julia, and GNU bash.

We have created a set of `Jupyter notebooks <http://jupyter.org/>`_
that you can run if you are interested in re-running partially or
totally our analysis.
In each notebook you will find further details for running them.

|

|

Requirements
************

To be able to run those notebooks you need to have the following
software installed:

- `Software requirements <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/software_requeriments.ipynb>`_: Instructions to install requeriments

|

|


How to reproduce the analyses
*****************************

We divided the analyses in different parts. Please look at the different jupyter notebooks which will point
to the commands needed for reproducing the analysis.

- `preprocessing <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/preprocessing_data.ipynb>`_: Creates and preprocesses some files.

- `formatting and annotating <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/formatting_and_annotating.ipynb>`_: Formats mutation files and annotates them.

- `signature extraction <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/signature_extraction.ipynb>`_: Extracts mutational signatures.

- `drugs preprocessing <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/drugs_preprocessing.ipynb>`_: Processing the drugs administered to the patients before the biopsy

- `regression analysis <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/regression.ipynb>`_: Perform the logistic regression on treatment-signatures.

- `mutational toxicity <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/mutational_toxicity.ipynb>`_: Calculate the mutational toxicity per sample.

- `timing analysis <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/timing_analysis.ipynb>`_: Get the early/late, clonal/subclonal activity per mutational signature.

- `leishmania treated analysis <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/leishmania_data_and_variant_calling.ipynb>`_: Analyze leishmania samples treated with 5-FU.

- `signature mixing <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/signature_mixing.ipynb>`_: Explore the putative components of mutational signatures.

- `extended_figures <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/extended_figures.ipynb>`_: Perform extra analyses presented in Extended Data Figures.


Each folder contains a notebook with a brief description and the requirements (notebooks that would need to be executed before).

If you would like to reproduce a specific figure, the following notebook will direct you to the analysis that generates it.

- `figures <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/figures.ipynb>`_: Locates each figure to the respective analysis notebook.

|

|

Simulations
***********

We also provide the code to reproduce the simulations performed in the supplementary materials `here <https://bitbucket.org/bbglab/mutfootprints/src/master/simulations/>`_.