# Import modules
import os
import sys
import gzip
import pickle

import pandas as pd
from tqdm import tqdm

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from utils.serializable import leveled_default_dict


# this will create a dictionary with the ML assignment for each signature
def ML_assignment_dictionary(type_mut, process_file):

    rev = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    name_tumor = os.path.basename(process_file).split('.processes')[0]
    outfile_name = '{}/{}.ML_assign.pckl.gz'.format(os.path.dirname(process_file), name_tumor)

    if not os.path.isfile(outfile_name):

        class_d = leveled_default_dict(str, level=2)
        class_full = leveled_default_dict(dict, level=2)

        exposure_file = process_file.replace('processes', 'exposures')

        if os.path.isfile(exposure_file):

            W = pd.read_csv(process_file, sep='\t', index_col=0)
            H = pd.read_csv(exposure_file, sep='\t', index_col=0)

            H = H.T.loc[:, ~H.T.columns.duplicated()].T
            W = W.loc[:, ~W.columns.duplicated()]

            for sample_index, sample in tqdm(enumerate(H.columns)):

                # multiply exposures by W
                ex = H[sample] * W
                ex['row_sum'] = ex.sum(axis=1)
                new = ex.div(ex['row_sum'], axis=0)
                new.drop('row_sum', axis=1, inplace=True)
                class_full[sample] = new.to_dict(orient='index')

                # get the one with the maximum probability
                res = new.idxmax(axis=1)

                for i, row in (W * H[sample]).iterrows():

                    # assign the signature to this context
                    class_d[sample][i] = res[i]

                    # now the complementary reverse too, for DBS and SNV
                    if type_mut == 'snvs':
                        rev_order = '{}[{}>{}]{}'.format(rev[i[-1]], rev[i[2]], rev[i[4]], rev[i[0]])
                        class_d[sample][rev_order] = res[i]

                    elif type_mut == 'dbs':
                        rev_order = '{}{}_{}{}'.format(rev[i[1]], rev[i[0]], rev[i[-1]], rev[i[-2]])
                        class_d[sample][rev_order] = res[i]

                for i, v in new.to_dict(orient='index').items():

                    if type_mut == 'snvs':
                        rev_order = '{}[{}>{}]{}'.format(rev[i[-1]], rev[i[2]], rev[i[4]], rev[i[0]])
                        class_full[sample][rev_order] = v

                    elif type_mut == 'dbs':
                        rev_order = '{}{}_{}{}'.format(rev[i[1]], rev[i[0]], rev[i[-1]], rev[i[-2]])
                        class_full[sample][rev_order] = v

            pickle.dump(dict(class_d), gzip.open(outfile_name, 'wb'))
            outfile_name = '{}/{}.{}.total_assign.pckl.gz'.format(os.path.dirname(process_file), name_tumor, type_mut)
            pickle.dump(dict(class_full), gzip.open(outfile_name, 'wb'))
        else:
            print(exposure_file, "not found")
    else:
        class_d = pickle.load(gzip.open(outfile_name))

    return class_d


# retrieve ML
def assing_ML(df, dML):

    return dML[df['SAMPLE']].get(df['VARIANT_CLASS'], 'not_found')


# provide
def doMLsignatures(path_original, path_processes_file, outpath):

    # create outpath if it doesnt exist
    os.makedirs(outpath, exist_ok=True)
    df = pd.read_csv(path_original, sep='\t', low_memory=False)

    name_ttype = os.path.basename(path_original).split('.')[0]
    type_mut = os.path.basename(path_original).split('.')[1]

    # get the dictionary with the MLikelihood
    dML = ML_assignment_dictionary(type_mut, path_processes_file)

    all_samples_extraction = set(list(dML.keys()))
    all_samples = set(df['SAMPLE'].tolist())

    # get only samples that were used for the extraction
    wanted_samples = [sample for sample in all_samples if sample in all_samples_extraction]
    df = df[df['SAMPLE'].isin(wanted_samples)]
    df['ML'] = df.apply(assing_ML, axis=1, args=(dML,))

    # this is the full file, we can also generate the file per signature to parallelize afterwards
    df.to_csv('{}/{}.ML.{}.gz'.format(outpath, name_ttype, type_mut), sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    path_original = sys.argv[1]
    path_processes_file = sys.argv[2]
    outpath = sys.argv[3]
    doMLsignatures(path_original, path_processes_file, outpath)
