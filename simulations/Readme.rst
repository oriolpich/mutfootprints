
README
======

This folder contains the code used for the simulation analysis of `PAPER <LINK TO THE PAPER>`_.
If you use this code in a publication, please cite:

.. admonition:: Citation
   :class: note

   Oriol Pich, Ferran Muiños, Abel Gonzalez-Perez, Nuria Lopez-Bigas, 
   The mutational footprints of cancer therapies, XX doi: `XX <LINK>`_


Running this software
*********************

These analysis have been performed using software in Python, R and GNU bash.

We have created a set of `Jupyter notebooks <http://jupyter.org/>`_
that you can run if you are interested in re-running partially or
totally the analysis.


Requirements
************

Keep in mind that to be able to run those notebooks you need to have the following
software installed (we also indicate the version so you can
reproduce the exact same results):

- `Software requirements <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/software_requeriments.ipynb>`_: Instructions to install requeriments

For the specific aim of running the simulation analyses you can also create a new
conda environment as follows:

.. code-block::

   $ conda env create -f simulations.yml

or embed the dependencies of **simulations.yml** in your running conda environment.

How to reproduce the analyses
*****************************

We divided the analyses in different parts. Please look at the different jupyter notebooks which will point
to the commands needed for reproducing the analysis. The details are enclosed in Supplementary Note 2.

- `1_synthetic_build <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/simulations/1_synthetic_build.ipynb>`_: Create synthetic mutational catalogues

- `2_deconstruction <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/2_deconstruction.ipynb>`_: Signature deconstruction and similarity against injected signal

- `3_drug_signature_dependency <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/3_drug_signature_dependency.ipynb>`_: Statistical dependencies between injected exposure and inferred signals

- `4_placebo_cotreatments <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/4_placebo_cotreatments.ipynb>`_: Disentangle true vs placebo-induced associations

- `5_exposure_recovery <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/mutfootprints/raw/master/5_exposure_recovery.ipynb>`_: Provides different measures of concordance between injected exposure and inferred exposures
