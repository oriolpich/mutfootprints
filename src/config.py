#!/usr/bin/env python


Conf = {
    "mappability_file": 'data/mappability/hg19.mappable.nochr.gz',
    "unmappability_file": 'data/mappability_files/hg19.notmappable.gz',
    "unmappability_file_hg38": 'data/mappability_files/hg38.notmappable.nochr.gz',

    "transcription_file": "data/asymmetry_files/strand_coordinates.bed.gz",
    "replication_file": 'data/asymmetry_files/replication_domain_20kb.bed.gz',

    "SBS_PCAWG": 'data/signatures_files/PCAWG/sigProfiler_SBS_signatures_2018_03_28.csv',
    "DBS_PCAWG": 'data/signatures_files/PCAWG/sigProfiler_DBS_signatures.csv',
    "ID_PCAWG": 'data/signatures_files/PCAWG/sigProfiler_ID_signatures.csv',

    "preTreatments": "data/clinical_data/preTreatments_fixed.tsv",
    "FDA": "data/clinical_data/product.txt",

    'treatment_specific_drug': 'data/clinical_data/hartwig_FDA_treatments_specific.pckl.gz',
    'treatment_FDA_drug': 'data/clinical_data/hartwig_FDA_treatments.pckl.gz',
    'treatment_specific_drug_stringent': 'data/clinical_data/hartwig_FDA_treatments_stringent.pckl.gz',

    # CHANGE TO YOUR OWN PATHS
    "path_metadata": 'YOUR_HMF_PATH/DR-024-update2_metadata.tsv',
    "age_metadata": 'YOUR_HMF_PATH/DR-024-update2_metadata.tsv',
    "path_predrugs_original": "YOUR_HMF_PATH/preDrugsByPatient.tsv",
    'Path_STJUDE':'YOUR_STJUDE_PATH/',
    'CGC_path': 'data/megabase_probability/YOUR_CGC_file' # we used --> Census_allMon Feb 18 17_56_21 2019.tsv
}
