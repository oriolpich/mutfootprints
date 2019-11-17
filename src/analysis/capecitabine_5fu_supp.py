# Import modules
import sys
import os
from collections import defaultdict
import gzip

import pandas as pd
import seaborn as sns
import dill as pickle
import matplotlib.pyplot as plt

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from utils.metadata import return_metadata
from utils.colors_ttype import return_colors
from utils.plots_utils import config_params
from config import Conf


def get_treated_not_treated():
    # get treate samples of each type
    treated_samples = pickle.load(gzip.open(Conf['treatment_specific_drug']))
    samples_FU = set()
    samples_capecitabine = set()

    for sample in treated_samples['FLUOROURACIL']['Pan']['YES']:
        if sample not in treated_samples['CAPECITABINE']['Pan']['YES']:
            samples_FU.add(sample)

    for sample in treated_samples['CAPECITABINE']['Pan']['YES']:
        if sample not in treated_samples['FLUOROURACIL']['Pan']['YES']:
            samples_capecitabine.add(sample)

    non_treated = treated_samples['FLUOROURACIL']['Pan']['NO'] & treated_samples['CAPECITABINE']['Pan']['NO']

    return non_treated, samples_FU, samples_capecitabine


def get_items_to_plot(samples_set, expos, sig):

    df = pd.read_csv(expos, sep='\t', index_col=0)
    df = df.T

    dic_primary_full, _ = return_metadata()
    colors_dict = return_colors()

    vals_to_plot = defaultdict(list)
    colors = defaultdict(list)
    count_exposed = defaultdict(int)
    count_total = defaultdict(int)

    for sample in samples_set:
        ttype = dic_primary_full[sample]
        if sample in df.index:
            sample_exp = df.loc[sample][sig]
            if sample_exp > 1:
                count_exposed[ttype] += 1
                vals_to_plot[ttype].append(sample_exp)
                colors[ttype].append(colors_dict[ttype])
            count_total[ttype] += 1

    return vals_to_plot, colors, count_exposed, count_total


def get_range_nums(l, v):
    return [v+(num / len(l)) for num in range(len(l))]


def plot_scatter_bars(expos, sig):

    non_treated, samples_FU, samples_capecitabine =get_treated_not_treated()
    vals_to_plot_FU, colors_FU, count_exposed_FU, count_total_FU = get_items_to_plot(samples_FU, expos, sig)
    vals_to_plot_capecitabine, colors_capecitabine, count_exposed_cape, count_total_cape = get_items_to_plot(
        samples_capecitabine, expos, sig
    )

    for ttype in ['Colon-Rectum', 'Breast']:
        config_params(font_size=7)

        fig, ax = plt.subplots(1, 1, figsize=(0.5, 1.25))
        plt.yscale('log')
        ax.set_ylabel('SBS Capecitabine/5-FU')

        range_nums = get_range_nums(vals_to_plot_FU[ttype], 0)
        ax.scatter(range_nums,
                   sorted(vals_to_plot_FU[ttype]), color=colors_FU[ttype], s=1, alpha=0.75)

        range_nums = get_range_nums(vals_to_plot_capecitabine[ttype], 1)
        ax.scatter(
            range_nums,
            sorted(vals_to_plot_capecitabine[ttype]),
            color=colors_capecitabine[ttype],
            s=1,
            alpha=0.75
        )

        print(len(vals_to_plot_FU[ttype]), len(vals_to_plot_capecitabine[ttype]), ttype)
        plt.xticks([0.5, 1.5], ['5-FU', 'Capecitabine'], rotation=90)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.ylim(0, 10000)
        plt.savefig('figures/sigmoid_treat_{}.svg'.format(ttype))
        plt.show()
        fig, ax = plt.subplots(1, 1, figsize=(0.5, 0.4))

        ax.bar(0, count_exposed_FU[ttype]/count_total_FU[ttype], color='#2c89a0ff')
        ax.bar(1, count_exposed_cape[ttype]/count_total_cape[ttype], color='#2c89a0ff')
        plt.ylim(0, 1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylabel('Proportion  of samples\nwith activity')
        plt.savefig('figures/bar_treat_{}.svg'.format(ttype))
        plt.show()


def second_plot(path_to_exposures, sig_cape, sig_17):

    df = pd.read_csv(path_to_exposures, sep='\t', index_col=0).T
    treated_dic = pickle.load(gzip.open(Conf['treatment_specific_drug']))

    treated_col = treated_dic['5-FU_CAPE']['Colon-Rectum']['YES']
    treated_br = treated_dic['5-FU_CAPE']['Breast']['YES']

    untreated_col = treated_dic['5-FU_CAPE']['Colon-Rectum']['NO']
    untreated_br = treated_dic['5-FU_CAPE']['Breast']['NO']

    exp_17_col = df[sig_17].loc[treated_col]
    exp_17_br = df[sig_17].loc[treated_br]

    exp_cape_col = df[sig_cape].loc[treated_col]
    exp_cape_br = df[sig_cape].loc[treated_br]

    notexp_17_col = df[sig_17].loc[untreated_col]
    notexp_17_br = df[sig_17].loc[untreated_br]

    notexp_cape_col = df[sig_cape].loc[untreated_col]
    notexp_cape_br = df[sig_cape].loc[untreated_br]

    fig, ax = plt.subplots(1, 1, figsize=(1.5, 1))
    sns.boxplot(
        data=[
            notexp_17_col, notexp_cape_col, exp_17_col, exp_cape_col,
            notexp_17_br, notexp_cape_br, exp_17_br, exp_cape_br,
        ],
        palette=[
            '#053061', '#053061', '#053061', '#053061',
            '#fa9fb5', '#fa9fb5', '#fa9fb5', '#fa9fb5',
        ],
        showfliers=False,
        boxprops=dict(alpha=.9, ),
        linewidth=0.6
    )

    print(
        len(notexp_17_col), len(notexp_cape_col),
        len(exp_17_col), len(exp_cape_col),
        len(notexp_17_br), len(notexp_cape_br),
        len(exp_17_br), len(exp_cape_br)
    )

    plt.ylim(0, 5000)
    plt.savefig('figures/comparision_17.svg')
    plt.show()


def run():

    path_siganalyzer = "data/hartwig/signatures/extraction/results/SignatureAnalyzer/"
    path_siganalyzer += "snvs/exposures/Pan_full/Pan_full.exposures.tsv"
    sig = '31_SBS17b_0.968799_1'
    sig_other_17 = '19_SBS17b_0.932022_1'

    plot_scatter_bars(path_siganalyzer, sig)
    second_plot(path_siganalyzer, sig, sig_other_17)


if __name__ == '__main__':
    run()
