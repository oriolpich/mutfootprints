import sys, os
scripts_path = os.path.abspath(os.path.join(__file__, '..', '..' ))
sys.path.insert(0, scripts_path)

from utils.metadata import return_metadata
from config import Conf
import dill as pickle
import gzip
import pandas as pd
from functools import partial
from contextlib import suppress

def create_samples_dict(sigprofiler_sigs, treatments, mutations ):

    samples_dict = {}
    for treatm in sigprofiler_sigs:
        if treatm not in ['AGING', 'SMOKING']:
            samples = set()
            for ttype in treatments[treatm]:
                samples = samples.union(treatments[treatm][ttype]['YES'])
            samples_dict[treatm] = samples
    samples_dict['AGING'] = set(mutations['sample'].values)
    samples_dict['SMOKING'] = mutations[mutations['ttype'] == 'Lung']['sample'].unique()

    with gzip.open('data/megabase_probability/samples_per_treatment.pickle.gz', 'wb') as f:
        pickle.dump(samples_dict, f)

    return samples_dict

def get_assignment(row, treatm, samples):

    if (row['sample'] in samples) and (row['sample'] in assign):
        return assign[row['sample']][row['context']][sigprofiler_sigs[treatm]]
    return 0

def create_probability_in_MB():

    mutations_path = 'data/megabase_probability/pan.snvs.no_cgc_mappable_nottype.gz'
    mutations = pd.read_csv(mutations_path, sep='\t', compression='gzip', header=None,
                            usecols=[5, 6, 16])
    mutations.columns = ['context', 'sample', 'megabase']
    mutations['context'] = mutations['context'].astype('category')
    mutations['sample'] = mutations['sample'].astype('category')
    mutations['megabase'] = mutations['megabase'].astype('category')

    treatments_path = Conf['treatment_specific_drug']
    with gzip.open(treatments_path, 'rb') as f:
        treatments = pickle.load(f)

    sample2ttype, _ = return_metadata()
    mutations['ttype'] = mutations['sample'].apply(lambda x: sample2ttype[x])
    mutations['ttype'] = mutations['ttype'].astype('category')

    samples_dict = create_samples_dict(sigprofiler_sigs, treatments, mutations )
    create_file_mutations_probability(samples_dict, mutations)
    create_samples_per_ttm_ttype(samples_dict)
    raw_megabase_mutations_ttype()
    normalized_prob()


def create_file_mutations_probability(samples_dict, mutations):

    path ="data/megabase_probability/"
    for treatm in ['CAPECITABINE', 'CARBOPLATIN', 'CISPLATIN', 'OXALIPLATIN']:
        samples = samples_dict[treatm]
        task = partial(get_assignment, treatm=treatm, samples=samples)
        mutations['probability'] = mutations.apply(task, axis=1)
        mutations.to_csv(path + treatm + '_mutations_' + 'probability.tsv.gz', sep='\t', index=False, compression='gzip')

    for treatm in ['AGING']:
        samples = set(mutations['sample'].values)
        task = partial(get_assignment, treatm=treatm, samples=samples)
        mutations['probability'] = mutations.apply(task, axis=1)
        mutations.to_csv(path + treatm + '_mutations_' + 'probability.tsv.gz', sep='\t', index=False, compression='gzip')

    for treatm in ['SMOKING']:
        samples = mutations[mutations['ttype'] == 'Lung']['sample'].unique()
        task = partial(get_assignment, treatm=treatm, samples=samples)
        mutations['probability'] = mutations.apply(task, axis=1)
        mutations.to_csv(path + treatm + '_mutations_' + 'probability.tsv.gz', sep='\t', index=False, compression='gzip')

def create_samples_per_ttm_ttype(samples_dict):

    path ="data/megabase_probability/"
    d = {}
    for treatm in ['CAPECITABINE', 'CARBOPLATIN', 'CISPLATIN', 'OXALIPLATIN', 'AGING', 'SMOKING']:
        d[treatm] = {}
        df = pd.read_csv(path + treatm + '_mutations_' + 'probability.tsv.gz', sep='\t', compression='gzip')
        g_per_ttype = df.groupby('ttype')
        for k, v in g_per_ttype:
            samples = set(list(v['sample'].unique())).intersection(samples_dict[treatm])
            d[treatm][k] = len(samples)

    with gzip.open(path +'samples_per_treatment_ttype.pickle.gz', 'wb') as f:
        pickle.dump(d, f)

def raw_megabase_mutations_ttype():

    path ="data/megabase_probability/"
    d = {}
    for treatm in ['CAPECITABINE', 'CARBOPLATIN', 'CISPLATIN', 'OXALIPLATIN', 'AGING', 'SMOKING']:
        d[treatm] = {}
        df = pd.read_csv(path + treatm + '_mutations_' + 'probability.tsv.gz', sep='\t', compression='gzip')
        g_per_ttype = df.groupby('ttype')
        for k, v in g_per_ttype:
            aggregated = v.groupby('megabase')['probability'].agg('sum')
            d[treatm][k] = dict(zip(aggregated.index, aggregated.values))
    with gzip.open(path + 'raw_megabase_mutations_per_ttype.pickle.gz', 'wb') as f:
        pickle.dump(d, f)

def normalized_prob():

    path ="data/megabase_probability/"

    raw_megabase = pickle.load(gzip.open(path +'raw_megabase_mutations_per_ttype.pickle.gz', 'rb'))
    samples_per_treatment = pickle.load(gzip.open(path + 'samples_per_treatment_ttype.pickle.gz', 'rb'))
    megabase_mappable = pickle.load(gzip.open(path + 'mappable_counts_megabase_mutations.pckl.gz', 'rb'))

    normalized_raw_megabase = {}
    for treatm in raw_megabase:

        normalized_raw_megabase[treatm] = {}
        for ttype in raw_megabase[treatm]:
            normalized_raw_megabase[treatm][ttype] = {}
            n_samples = samples_per_treatment[treatm][ttype]
            for chunk in raw_megabase[treatm][ttype]:
                print(chunk, n_samples, treatm)
                k = n_samples * megabase_mappable[chunk]
                if k > 0:
                    normalized_raw_megabase[treatm][ttype][chunk] = raw_megabase[treatm][ttype][chunk] / k
                else:
                    normalized_raw_megabase[treatm][ttype][chunk] = 0

    with gzip.open(path + 'normalized_megabase_mutations_ttype.pickle.gz', 'wb') as f:
        pickle.dump(normalized_raw_megabase, f)


assign_path = 'data/hartwig/signatures/extraction/results/SigProfiler/snvs/processes/PanNoSkinNoUnknown.snvs/PanNoSkinNoUnknown.snvs.snvs.total_assign.pckl.gz'
with gzip.open(assign_path, 'rb') as f:
    assign = pickle.load(f)

sigprofiler_sigs = {'CAPECITABINE': '19_SBS17b_0.961548_0.99',
                    'CARBOPLATIN': '1_SBS31_0.968153_0.98',
                    'CISPLATIN': '1_SBS31_0.968153_0.98',
                    'OXALIPLATIN': '20_0.92',
                    'AGING': '23_SBS1_0.996494_1.0',
                    'SMOKING': '17_SBS4_0.935484_0.97'}

create_probability_in_MB()

