# Import modules
import os
import sys
from collections import defaultdict
import gzip
from multiprocessing import Pool

import pandas as pd
from tqdm import tqdm
import dill as pickle

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf
from select_signatures import ensemble_regression


# Global variable
NCPUS = 4  # number of CPUs for multiprocessing


def do_regression(s):

    drug, type_mut, tumor_type, outpath, treated_samples_file, path_to_exposures = s

    try:
        effect, pvals, merge_res, treated_samples = ensemble_regression(
            drug, treated_samples_file, path_to_exposures, ttype=tumor_type, n=1000, keep_sigs=False
        )

        # get the results and save them
        dfe = defaultdict(dict)
        dfe['effect_size'] = effect
        dfe['pvals'] = pvals
        dfe['treatment'] = str.lower(drug.replace(' ', '_').replace('/', '_'))

        dfe['ttype'] = tumor_type

        out = pd.DataFrame(dfe)
        out['signature'] = out.index.tolist()
        out['total_treated'] = treated_samples

        # prepare the outputs
        outpath_mut = '{}/{}'.format(outpath, type_mut)
        os.makedirs(outpath_mut, exist_ok=True)
        type_regression = os.path.basename(treated_samples_file).split('.')[0]
        formatted_drug = drug.replace(' ', '_').replace('/', '_').replace('-', '_').replace('(', '_').replace(')', '_')

        outfile = '{}/{}.{}.{}.tsv'.format(outpath_mut, tumor_type, formatted_drug,  type_regression)
        out.to_csv(outfile, sep='\t', index=False, header=True)
    except Exception:
        pass


# get the subset of tumor type - drug that will be worth to explore
def get_drugs_prevalent_specific(treated_file):

    # load specifically treated samples
    treated_samples = pickle.load(gzip.open(treated_file))

    # get the amount of patients treated per drug/ttype
    countTreated = defaultdict(lambda: defaultdict(int))
    for drug, data in treated_samples.items():
        for ttype, dic in data.items():
            if ttype != 'Pan':
                treated_in_ttype = len(dic['YES'])
                countTreated[drug][str(ttype)] = treated_in_ttype

    total_count_df = pd.DataFrame(countTreated)
    norm = (total_count_df / total_count_df.sum()).fillna(0)
    assoc = set()
    prevalent = set()

    # select those cases where more than the 25% of the treated patients come from one tumor type
    for col in norm.columns:
        for ttype, row in norm.iterrows():
            if row[col] > 0.25:
                assoc.add((col, ttype))
                prevalent.add(col)

    if len(total_count_df) > 1:

        # the ttype with more patients treated with each drug
        max_patients = total_count_df.idxmax().to_dict()
        set_enough_samples = set()
        for drug, ttype in max_patients.items():

            val = total_count_df[drug].loc[ttype]
            if (val > 25) & (drug in prevalent):
                set_enough_samples.add((drug, ttype))
    else:
        set_enough_samples = {}
        max_patients = {}

    return prevalent, set_enough_samples, max_patients


def multiprocess_regression(treated_file, path_exposures, outpath, type_mut, ncpu):

    treated_dict = pickle.load(gzip.open(treated_file))
    prevalent, set_enough_samples, max_patients = get_drugs_prevalent_specific(treated_file)

    # get the type of mutation

    # PAN ANALYSIS
    tumor_type = "Pan"
    all_drugs = [
        (drug, type_mut, tumor_type, outpath, treated_file, path_exposures)
        for drug in treated_dict.keys() if len(treated_dict[drug][tumor_type]['YES']) > 25
    ]

    with Pool(ncpu) as pool:
        for _ in tqdm(pool.imap_unordered(do_regression, all_drugs), total=len(all_drugs)):
            continue

    # TUMOR TYPE PREVALENT ANALYSIS
    if len(set_enough_samples) > 1:
        all_drugs = [
            (drug, type_mut, tumor_type, outpath, treated_file, path_exposures)
            for drug, tumor_type in set_enough_samples
        ]

        with Pool(ncpu) as pool:
            for _ in tqdm(pool.imap_unordered(do_regression, all_drugs), total=len(all_drugs)):
                continue


def run():

    treated_file_specific_drug = Conf['treatment_specific_drug']
    treated_file_FDA = Conf['treatment_FDA_drug']
    todo = [treated_file_specific_drug, treated_file_FDA]

    outpath = 'data/regression_treatments/SignatureAnalyzer/'
    os.makedirs(outpath, exist_ok=True)

    for treated_path in todo:

        # --------
        # SNVS
        # --------
        path_exposures = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/snvs/exposures/Pan_full/Pan_full.exposures.tsv'
        type_mut = 'snv'
        multiprocess_regression(treated_path, path_exposures, outpath, type_mut,  NCPUS)

        # -------
        # DBS
        # -------
        path_exposures = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/dbs/exposures/Pan/Pan.exposures.tsv'
        type_mut = 'dbs'
        multiprocess_regression(treated_path, path_exposures, outpath, type_mut, NCPUS)

        # --------
        # INDELS
        # --------
        path_exposures = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/indels/exposures/Pan/Pan.exposures.tsv'
        type_mut = 'indels'
        multiprocess_regression(treated_path, path_exposures, outpath, type_mut,  NCPUS)

    outpath = 'data/regression_treatments/SigProfiler/'
    os.makedirs(outpath, exist_ok=True)
    base_path = 'data/hartwig/signatures/extraction/results/SigProfiler/'

    for treated_path in todo:
        # --------
        # SNVS
        # --------
        path_exposures = base_path + 'snvs/exposures/PanNoSkinNoUnknown.snvs/PanNoSkinNoUnknown.snvs.exposures.tsv'
        type_mut = 'snv'
        multiprocess_regression(treated_path, path_exposures, outpath, type_mut,  NCPUS)

        # -------
        # DBS
        # -------
        path_exposures = base_path + 'dbs/exposures/PanNoSkinNoUnknown.dbs/PanNoSkinNoUnknown.dbs.exposures.tsv'
        type_mut = 'dbs'
        multiprocess_regression(treated_path, path_exposures, outpath, type_mut, NCPUS)

        # --------
        # INDELS
        # --------
        path_exposures = base_path + 'indels/exposures/PanNoSkinNoUnknown.indels/PanNoSkinNoUnknown.indels.exposures.tsv'
        type_mut = 'indels'
        multiprocess_regression(treated_path, path_exposures, outpath, type_mut,  NCPUS)


if __name__ == '__main__':
    run()
