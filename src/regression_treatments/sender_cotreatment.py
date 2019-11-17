# Import modules
import os
import sys
import pickle
import gzip
from collections import defaultdict
from glob import glob

import pandas as pd

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf
from cotreatment import attribution


def merge_specific(folder, outpath):
    os.makedirs(outpath, exist_ok=True)
    for mut in ['snv', 'dbs', 'indels']:
        pool = []
        for fn in glob(os.path.join(folder, mut, "*.hartwig_FDA_treatments_specific.tsv")):
            df = pd.read_csv(fn, sep='\t')
            pool.append(df)
        if len(pool):
            df = pd.concat(pool)
            df.to_csv(
                os.path.join(outpath, 'regression_results_treatments_specific_{}.tsv'.format(mut)),
                sep='\t',
                index=False
            )
            pool = []
            for fn in glob(os.path.join(folder, mut, "*.hartwig_FDA_treatments.tsv")):
                df = pd.read_csv(fn, sep='\t')
                pool.append(df)
            df = pd.concat(pool)
            df.to_csv(
                os.path.join(outpath, 'regression_results_treatments_FDA_{}.tsv'.format(mut)),
                sep='\t',
                index=False
            )


def create_matrix_treatment(treated_samples_specific):

    treated_samples = pickle.load(gzip.open(treated_samples_specific))

    d_fill = defaultdict(dict)
    for drug in treated_samples:
        for sample in treated_samples[drug]['Pan']['YES']:
            formatted_drug = str.lower(drug.replace(' ', '_').replace('/', '_'))
            d_fill[formatted_drug][sample] = 1
    new_matrix_treatment = pd.DataFrame(d_fill).fillna(0)

    path_file = 'data/clinical_data/matrix_{}.tsv'.format(os.path.basename(treated_samples_specific).split('.')[0])
    new_matrix_treatment.to_csv(
        path_file, sep='\t', index=True, header=True
    )

    return path_file


def run():
    # ---------------
    # SigProfiler
    # ---------------
    folder = 'data/regression_treatments/SigProfiler/'
    outpath = 'data/regression_treatments/SigProfiler/merged/'
    merge_specific(folder, outpath)
    path_exposures_base = 'data/hartwig/signatures/extraction/results/SigProfiler/{}/exposures/'
    path_exposures_base += 'PanNoSkinNoUnknown.{}/PanNoSkinNoUnknown.{}.exposures.tsv'

    # SNVS
    matrix_treatment_path_full = create_matrix_treatment(Conf['treatment_FDA_drug'])
    matrix_treatment_path_specific = create_matrix_treatment(Conf['treatment_specific_drug'])

    path_exposures = path_exposures_base.format('snvs', 'snvs', 'snvs')
    discriminate_path_full = outpath + 'regression_results_treatments_FDA_snv.tsv'
    discriminate_path_specific = outpath + 'regression_results_treatments_specific_snv.tsv'

    output_full = outpath + 'cotreatment_FDA_snvs.tsv'
    output_specific = outpath + 'cotreatment_specific_drug_snvs.tsv'

    attribution(path_exposures, discriminate_path_full, matrix_treatment_path_full, output_full)
    attribution(path_exposures, discriminate_path_specific, matrix_treatment_path_specific, output_specific)

    # DBS
    path_exposures = path_exposures_base.format('dbs', 'dbs', 'dbs')
    discriminate_path_full = outpath + 'regression_results_treatments_FDA_dbs.tsv'
    discriminate_path_specific = outpath + 'regression_results_treatments_specific_dbs.tsv'

    output_full = outpath + 'cotreatment_FDA_dbs.tsv'
    output_specific = outpath + 'cotreatment_specific_drug_dbs.tsv'

    attribution(path_exposures, discriminate_path_full, matrix_treatment_path_full, output_full)
    attribution(path_exposures, discriminate_path_specific, matrix_treatment_path_specific, output_specific)

    # INDELS
    path_exposures = path_exposures_base.format('indels', 'indels', 'indels')
    discriminate_path_full = outpath + 'regression_results_treatments_FDA_indels.tsv'
    discriminate_path_specific = outpath + 'regression_results_treatments_specific_indels.tsv'

    output_full = outpath + 'cotreatment_FDA_indels.tsv'
    output_specific = outpath + 'cotreatment_specific_drug_indels.tsv'

    attribution(path_exposures, discriminate_path_full, matrix_treatment_path_full, output_full)
    attribution(path_exposures, discriminate_path_specific, matrix_treatment_path_specific, output_specific)

    # ---------------
    # SignatureAnalyzer
    # ---------------
    folder = 'data/regression_treatments/SignatureAnalyzer/'
    outpath = 'data/regression_treatments/SignatureAnalyzer/merged/'
    merge_specific(folder, outpath)
    path_exposures_base = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/'

    # SNVS
    matrix_treatment_path_full = create_matrix_treatment(Conf['treatment_FDA_drug'])
    matrix_treatment_path_specific = create_matrix_treatment(Conf['treatment_specific_drug'])

    path_exposures = path_exposures_base + 'snvs/exposures/Pan_full/Pan_full.exposures.tsv'
    discriminate_path_full = outpath + 'regression_results_treatments_FDA_snv.tsv'
    discriminate_path_specific = outpath + 'regression_results_treatments_specific_snv.tsv'

    output_full = outpath + 'cotreatment_FDA_snvs.tsv'
    output_specific = outpath + 'cotreatment_specific_drug_snvs.tsv'

    attribution(path_exposures, discriminate_path_full, matrix_treatment_path_full, output_full)
    attribution(path_exposures, discriminate_path_specific, matrix_treatment_path_specific, output_specific)

    # DBS
    path_exposures = path_exposures_base + 'dbs/exposures/Pan/Pan.exposures.tsv'
    discriminate_path_full = outpath + 'regression_results_treatments_FDA_dbs.tsv'
    discriminate_path_specific = outpath + 'regression_results_treatments_specific_dbs.tsv'

    output_full = outpath + 'cotreatment_FDA_dbs.tsv'
    output_specific = outpath + 'cotreatment_specific_drug_dbs.tsv'

    attribution(path_exposures, discriminate_path_full, matrix_treatment_path_full, output_full)
    attribution(path_exposures, discriminate_path_specific, matrix_treatment_path_specific, output_specific)

    # INDELS
    path_exposures = path_exposures_base + 'indels/exposures/Pan/Pan.exposures.tsv'
    discriminate_path_full = outpath + 'regression_results_treatments_FDA_indels.tsv'
    discriminate_path_specific = outpath + 'regression_results_treatments_specific_indels.tsv'

    output_full = outpath + 'cotreatment_FDA_indels.tsv'
    output_specific = outpath + 'cotreatment_specific_drug_indels.tsv'

    attribution(path_exposures, discriminate_path_full, matrix_treatment_path_full, output_full)
    attribution(path_exposures, discriminate_path_specific, matrix_treatment_path_specific, output_specific)


if __name__ == '__main__':
    run()
