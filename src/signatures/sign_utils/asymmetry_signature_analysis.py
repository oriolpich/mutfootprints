# Import modules
import sys
import os
from collections import defaultdict
from glob import glob

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.spatial.distance as spdist

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..', '..'))
sys.path.insert(0, scripts_path)


from utils.plots_utils import config_params


# split into even chunks
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


def get_cosinus_similarity_extractions(repli, trans):

    df_repli = pd.read_csv(repli, sep='\t', index_col=0)
    df_trans = pd.read_csv(trans, sep='\t', index_col=0)

    dic_trans = defaultdict(list)
    for signature in df_trans:
        val_sum = []
        sig_full = df_trans[signature].tolist()
        for ix, c in enumerate(chunks(sig_full, 2)):
            val_sum.append(np.sum(c))
        dic_trans[signature] = val_sum

    dic_repli = defaultdict(list)
    for signature in df_repli:
        val_sum = []
        sig_full = df_repli[signature].tolist()
        for ix, c in enumerate(chunks(sig_full, 2)):
            val_sum.append(np.sum(c))
        dic_repli[signature] = val_sum

    # find similar signatures to those reported in SigProfiler
    cos_sim = defaultdict(dict)
    for sig, vec1 in dic_repli.items():
        for s, vec2 in dic_trans.items():
            c = 1 - round(spdist.cosine(vec1, vec2), 7)
            cos_sim[sig][s] = c

    # select those with higher similarity per each of the signatures
    cos_df = pd.DataFrame(cos_sim)
    index_max = cos_df.idxmax()
    vals = cos_df.max()
    similar = pd.DataFrame(list(zip(index_max, vals)))

    similar.index = cos_df.columns
    high_cos = similar[similar[1] > 0.95]
    dic_equi = dict(zip(high_cos.index, high_cos[0]))

    print(dic_equi)
    return dic_equi


def return_values_components(all_f, dic_equi):

    axis_transcription = []
    axis_replication = []
    selected_sigs = defaultdict(list)
    for f in all_f:
        sig_rep = (os.path.basename(f).split('replication.')[1].split('.max_com')[0])
        if sig_rep in dic_equi:
            print(sig_rep)
            new_f = f.replace('replication', 'transcription').replace(sig_rep, dic_equi[sig_rep])
            df_rep = pd.read_csv(f, header=None, sep='\t')
            val_rep = 2*abs(0.5 - float(df_rep[1]/(df_rep[1]+df_rep[2])))
            df_trans = pd.read_csv(new_f, header=None, sep='\t')
            val_trans = 2*abs(0.5 - float(df_trans[1]/(df_trans[1] + df_trans[2])))
            axis_transcription.append(val_trans)
            axis_replication.append(val_rep)
            selected_sigs[sig_rep] = (val_rep, val_trans,)

    return axis_replication, axis_transcription, selected_sigs


def do_plot(type_extraction, axis_replication, axis_transcription, selected_sigs, wanted_sigs):

    config_params(6.5)
    fig, ax = plt.subplots(1, 1, figsize=(1.85, 2.1,))
    plt.scatter(axis_replication, axis_transcription, s=1, c='grey')

    for sig in wanted_sigs:
        print(sig, selected_sigs[sig])
        x, y = selected_sigs[sig]

        plt.scatter(x, y, s=10, c='#ff6600ff')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.ylabel('Transcription strand asymmetry')
    plt.xlabel('Replication strand asymmetry')

    plt.savefig('figures/{}/repli_trans_assym.svg'.format(type_extraction))
    plt.show()


def run():
    repli = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/snvs/processes/Pan_replication/Pan_replication.processes.tsv'
    trans = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/snvs/processes/Pan_transcription/Pan_transcription.processes.tsv'
    wanted = [
        '45_SBS7a_0.911948_1', '19_SBS13_0.942988_1', '17_SBS17b_0.961721_1',
        '4_SBS1_0.994881_1', '24_SBS4_0.965533_1', '21_SBS31_0.929973_1', '3_1'
    ]

    # get signature equivalences
    dic_equi = get_cosinus_similarity_extractions(repli, trans)

    all_f = glob('data/hartwig/signatures/extraction/results/SignatureAnalyzer/snvs/processes/Pan_replication/*max_*')
    axis_replication, axis_transcription, selected_sigs = return_values_components(all_f, dic_equi)

    type_extraction = 'SignatureAnalyzer'
    do_plot(type_extraction, axis_replication, axis_transcription, selected_sigs, wanted)

    repli = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/dbs/processes/Pan_replication/Pan_replication.processes.tsv'
    trans = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/dbs/processes/Pan_transcription/Pan_transcription.processes.tsv'
    wanted = ['4_DBS5_0.928311_1']

    # get signature equivalences
    dic_equi = get_cosinus_similarity_extractions(repli, trans)

    all_f = glob('data/hartwig/signatures/extraction/results/SignatureAnalyzer/dbs/processes/Pan_replication/*max_*')
    axis_replication, axis_transcription, selected_sigs = return_values_components(all_f, dic_equi)

    type_extraction = 'SignatureAnalyzer'
    do_plot(type_extraction, axis_replication, axis_transcription, selected_sigs, wanted)


if __name__ == '__main__':
    run()
