""" Load and preprocess signatures and signature activities instances """

import pandas as pd

import utils


# config paths

sa_signatures_path = './data/SignatureAnalyzer_COMPOSITE_SBS_W96.signature.031918.txt'
sa_activities_path = './data/SignatureAnalyzer_COMPOSITE_SNV.activity.FULL_SET.031918.txt'
sp_signatures_path = './data/sigProfiler_SBS.csv'


def load_preprocess_sa():

    signatures_sa = pd.read_csv(sa_signatures_path, sep='\t')

    activities_sa = pd.read_csv(sa_activities_path, sep='\t', index_col=0)

    signatures_sa.rename(columns=dict(list(map(lambda x: (x, x.split('_')[-2]),
                                           list(activities_sa.index)))),
                         inplace=True)

    signatures_sa['feature'] = signatures_sa.apply(lambda r: utils.key_func2(r), axis=1)
    signatures_sa.set_index('feature', inplace=True)

    signatures_sa_dict = {sig: dict(zip(list(signatures_sa.index),
                               list(signatures_sa[sig].values))) for sig in signatures_sa.columns}

    activities_sa.rename(index=dict(list(map(lambda x: (x, x.split('_')[-2]),
                                    list(activities_sa.index)))),
                         inplace=True)

    return signatures_sa, activities_sa, signatures_sa_dict


def load_preprocess_sp():

    signatures_sp = pd.read_csv(sp_signatures_path, sep=',')

    signatures_sp['key'] = signatures_sp.apply(lambda r: utils.key_func1(r), axis=1)
    signatures_sp.set_index('key', inplace=True)

    del signatures_sp['Type']
    del signatures_sp['SubType']

    signatures_sp_dict = {sig: dict(zip(list(signatures_sp.index),
                               list(signatures_sp[sig].values))) for sig in signatures_sp.columns}

    return signatures_sp, signatures_sp_dict


def retrieve_activities(tumor_types, activities, signatures):

    pool_cols = list(filter(lambda x: x.split('__')[0] in tumor_types, activities.columns))
    pool = activities.reindex(index=signatures.columns, fill_value=0)
    pool = pool.loc[signatures.columns, pool_cols]
    pool.fillna(0, inplace=True)

    return pool
