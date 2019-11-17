# Improt modules
import sys
import os
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..', '..'))
sys.path.insert(0, scripts_path)

from utils.metadata import return_metadata
from utils.plots_utils import config_params


def read_exposures(exposures_path):
    df = pd.read_csv(exposures_path, sep='\t', index_col=0)
    df = df.T
    return df


def plot_piecharts_signatures(exposures_path, type_mut,  type_extraction, figsize, min_val):

    colors_full = {
        "SBS": {
            25: '#003d7c',
            10: '#224ba5',
            5: '#4459ce',
            2.5: '#6767f7',
            1: '#9968c8',
            0.5: '#cc6999',
            0.25: '#ff6b6b',
            0.1: '#ff8962',
            0.05: '#ffa759',
            0: '#ffc651'
        },
        "ID": {
            1: '#003d7c',
            0.5: '#224ba5',
            0.1: '#4459ce',
            0.05: '#6767f7',
            0.04: '#9968c8',
            0.03: '#cc6999',
            0.02: '#ff6b6b',
            0.01: '#ff8962',
            0.001: '#ffa759',
            0: '#ffc651'
        },
        "DBS": {
            0.5: '#003d7c',
            0.1: '#224ba5',
            0.05: '#4459ce',
            0.03: '#6767f7',
            0.02: '#9968c8',
            0.01: '#cc6999',
            0.008: '#ff6b6b',
            0.005: '#ff8962',
            0.001: '#ffa759',
            0: '#ffc651'
        }
    }

    colors = colors_full[type_mut]
    result = read_exposures(exposures_path)
    config_params()
    dic_primary_full, _ = return_metadata()

    # result = pd.concat([df, df_mela])
    result = result.fillna(0)
    signatures = result.columns.tolist()

    # get list of similar found signatures in the extraction
    similar = [s for s in signatures if type_mut in s]
    notsimilar = [s for s in signatures if type_mut not in s]

    result['TTYPE'] = [dic_primary_full[t] for t in result.index.tolist()]

    dic_sig = defaultdict(lambda: defaultdict(float))
    dic_proportion = defaultdict(lambda: defaultdict(float))
    for ttype, data in result.groupby(by='TTYPE'):
        data2 = data.copy()
        data2.drop('TTYPE', axis=1, inplace=True)
        for col in data2:
            # we normalize it by the number of MB in the human genome (3234)
            dic_sig[ttype][col] = data2[data2[col] > min_val][col].median() / 3234

            if type_mut not in col:
                dic_proportion[ttype][col.split('_')[0]] = len(data2[data2[col] > min_val]) / len(data2)
            else:
                dic_proportion[ttype][col] = len(data2[data2[col] > min_val]) / len(data2)

    medians = pd.DataFrame.from_dict(dic_sig)

    # sorting the already known signatures
    keep_order_similar = defaultdict(list)
    for s in similar:
        number = s.split('_')[1].split(type_mut)[1]
        try:
            keep_n = int(number)
        except Exception:
            keep_n = int(number[:-1])

        keep_order_similar[keep_n].append(str(s))

    order_prev_labels = []
    order_prev = []
    for i in sorted(keep_order_similar, reverse=True):
        all_s = []
        d_equiv = defaultdict(str)
        for sig in keep_order_similar[i]:
            ID_sig = '{}_{}_{}'.format(sig.split('_')[1], sig.split('_')[2], sig.split('_')[0])
            d_equiv[ID_sig] = sig
        sorted_final_k = sorted(d_equiv)
        sorted_sigs_list = [d_equiv[ss] for ss in sorted_final_k]

        for sim, sig in reversed(list(enumerate(sorted_sigs_list, start=1))):
            if len(keep_order_similar[i]) == 1:
                order_prev_labels.append('E-{}{} ({}-like, {})'.format(
                    type_mut,
                    sig.split('_')[0],
                    sig.split('_')[1],
                    round(float(sig.split('_')[2]), 3))
                )
                order_prev.append(sig)
            else:
                order_prev_labels.append('E-{}{} ({}-like {}, {})'.format(
                    type_mut,
                    sig.split('_')[0],
                    sig.split('_')[1],
                    sim,
                    round(float(sig.split('_')[2]), 3))
                )
                order_prev.append(sig)

    no_similar_signatures = medians.loc[notsimilar]
    new_index = [int(l.split('_')[0]) for l in no_similar_signatures.index.tolist()]
    no_similar_signatures.index = new_index
    no_similar_signatures.sort_index(inplace=True, ascending=True)

    names_notsimilar = ['E-{} {}'.format(type_mut, c) for c in no_similar_signatures.index.tolist()[::-1]]

    # merge new and old
    merged = pd.concat([no_similar_signatures.sort_index(ascending = False), medians.loc[order_prev], ])

    # merged = merged.loc[order_prev+small_newset.index.tolist()]
    merged_labels = names_notsimilar + order_prev_labels

    config_params(5)

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    # plt.grid(b=True, which='major',)

    for yval, (i, row) in enumerate(merged.iterrows()):
        for xval, t in enumerate(merged.columns.tolist()):
            val = row[t]
            if val > 0:
                color = None
                for number in sorted(colors.keys(), reverse=True):
                    if (val > number) & (color is None):
                        color_scatter = colors[number]
                        break
                if type_mut in str(i):
                    plt.scatter(xval, yval, c=color_scatter, s=dic_proportion[t][i] * 20)
                else:
                    plt.scatter(xval, yval, c=color_scatter, s=dic_proportion[t][str(i)] * 20)

    ax.set_xticks(np.arange(len(merged.T)))
    ax.set_xticklabels(merged.columns.tolist(), rotation=90)
    ax.set_yticks(np.arange(len(merged)))

    ax.set_yticklabels(merged_labels)
    ax.xaxis.set_ticks_position('top')
    ax.set_axisbelow(True)

    ax.yaxis.grid(color='gray', linestyle='dashed', alpha=0.3)
    ax.xaxis.grid(color='gray', linestyle='dashed', alpha=0.3)

    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    plt.ylim(-1, len(merged))
    plt.tight_layout()
    plt.savefig('figures/{}/supp1_{}.svg'.format(type_extraction, type_mut))
    plt.savefig('figures/{}/supp1_{}.png'.format(type_extraction, type_mut), dpi=600)

    plt.close()


def run():
    snvs_figsize = (4, 8.75)
    dbs_figsize = (3.75, 4.25)
    indels_figsize = (3.75, 4.25)

    # ---
    # SigProfiler
    # ---

    type_extraction = 'SigProfiler'

    min_val = 200
    exposures_path = 'data/hartwig/signatures/extraction/results/SigProfiler/snvs/exposures/PanNoSkinNoUnknown.snvs/PanNoSkinNoUnknown.snvs.exposures.tsv'
    plot_piecharts_signatures(exposures_path, 'SBS',  type_extraction, snvs_figsize, min_val)

    min_val = 20
    exposures_path = 'data/hartwig/signatures/extraction/results/SigProfiler/indels/exposures/PanNoSkinNoUnknown.indels/PanNoSkinNoUnknown.indels.exposures.tsv'
    plot_piecharts_signatures(exposures_path, "ID", type_extraction, indels_figsize, min_val)

    min_val = 20
    exposures_path = 'data/hartwig/signatures/extraction/results/SigProfiler/dbs/exposures/PanNoSkinNoUnknown.dbs/PanNoSkinNoUnknown.dbs.exposures.tsv'
    plot_piecharts_signatures(exposures_path, "DBS", type_extraction, dbs_figsize, min_val)

    # ---
    # SignatureAnalyzer
    # ---
    type_extraction = 'SignatureAnalyzer'

    min_val = 200
    exposures_path = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/snvs/exposures/Pan_full/Pan_full.exposures.tsv'
    plot_piecharts_signatures(exposures_path, 'SBS',  type_extraction, snvs_figsize, min_val)

    min_val = 20
    exposures_path = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/indels/exposures/Pan/Pan.exposures.tsv'
    plot_piecharts_signatures(exposures_path, "ID", type_extraction, indels_figsize, min_val)

    min_val = 20
    exposures_path = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/dbs/exposures/Pan/Pan.exposures.tsv'
    plot_piecharts_signatures(exposures_path, "DBS", type_extraction, indels_figsize, min_val)


if __name__ == '__main__':
    run()
