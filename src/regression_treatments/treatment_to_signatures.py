# Import modules
import os
import sys
import gzip
import pickle
from collections import defaultdict
from glob import glob

import pandas as pd

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)


def create_dics_association(method_extraction):

    dic_FDA_to_specific = pickle.load(gzip.open('data/clinical_data/hartwig_typeDrug.pckl.gz'))
    dic_FDA_to_specific_fix_names = defaultdict(list)

    for group, d in dic_FDA_to_specific.items():
        fixed_lower = ['_'.join(list(map(str.lower, drug.split(' ')))) for drug in d]
        dic_FDA_to_specific_fix_names['_'.join(group.lower().split(' '))] = fixed_lower

    allf = glob('data/regression_treatments/{}/merged/cotreatment_FDA*tsv'.format(method_extraction))
    toconcat = []

    for f in allf:
        df = pd.read_csv(f, sep='\t')
        df['type'] = f.split('_')[-1].split('.tsv')[0]
        df['real_sig'] = df['signature'] + "_" + df['type']
        toconcat.append(df)
    merged_FDA_group = pd.concat(toconcat)

    dic_group = defaultdict(set)
    for i, row in merged_FDA_group.iterrows():
        for specific_treatment in dic_FDA_to_specific_fix_names[row['origin']]:
            dic_group[row['real_sig']].add(specific_treatment)

    allf = glob('data/regression_treatments/{}/merged/cotreatment_specific*tsv'.format(method_extraction))
    toconcat = []
    final_dict = defaultdict(set)

    for f in allf:
        df = pd.read_csv(f, sep='\t')
        df['type'] = f.split('_')[-1].split('.tsv')[0]
        df['real_sig'] = df['signature'] + "_" + df['type']
        toconcat.append(df)
    merged_specific_group = pd.concat(toconcat)

    for signature, data in merged_specific_group.groupby(by='real_sig'):
        if signature in dic_group:
            for i, row in data.iterrows():
                if row['treatment'] in dic_group[signature]:
                    final_dict[row['treatment']].add(signature)
        else:
            for i, row in data.iterrows():
                if row['origin'] != 'None':
                    final_dict[row['origin']].add(signature)

    return final_dict


def run():

    drugs2sig = defaultdict(dict)
    for method in ['SigProfiler', 'SignatureAnalyzer']:
        drugs2sig[method] = create_dics_association(method)

    # we will select the cases where both drugs and the type of significant signatures are found in both methods to be more robust.
    print(drugs2sig['SigProfiler'])
    print(drugs2sig['SignatureAnalyzer'])


if __name__ == '__main__':
    run()
