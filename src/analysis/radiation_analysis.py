import os
import sys
import gzip


import dill as pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf
from utils.metadata import return_metadata
from utils.colors_ttype import return_colors
from utils.plots_utils import config_params


def load_exposures(f):
    d = pd.read_csv(f, sep='\t', index_col=0)

    return d.T


def sci_notation(number, sig_fig=2):
    ret_string = "{0:.{1:d}e}".format(number, sig_fig)

    a, b = ret_string.split("e")
    b = int(b)  # removed leading "+" and strips leading zeros too.
    string = a + r"x$10^{" + str(b) + "}$"

    return string


def select_HRdeficient_samples(exp_snv):

    # select samples with high Bracaness and low Bracaness
    low_braca = np.percentile(exp_snv['13_SBS3_0.870929_1'].tolist(), 25)
    high_braca = np.percentile(exp_snv['13_SBS3_0.870929_1'].tolist(), 75)
    samples_braca = exp_snv[exp_snv['13_SBS3_0.870929_1'] >= high_braca].index.tolist()
    samples_lowbraca = exp_snv[exp_snv['13_SBS3_0.870929_1'] <= low_braca].index.tolist()

    return samples_braca, samples_lowbraca


def select_treated_not_treated(treated_samples_full, treatment, breaks, exp):
    # remove treated samples with some compounds that will cause DBS too

    radiated_samples = [s for s in treated_samples_full[treatment]['Pan']['YES'] if s in exp.index]

    only_treated = set()
    for sample in radiated_samples:
        forbidden = False
        for breaking in breaks:
            if breaking != treatment:
                if sample in treated_samples_full[breaking]['Pan']['YES']:
                    forbidden = True
        if forbidden is not True:
            only_treated.add(sample)

    non_breakable_samples = [s for s in treated_samples_full[treatment]['Pan']['NO'] if s in exp.index]
    not_breaked = set()

    for sample in non_breakable_samples:

        forbidden = False
        for breaking in breaks:
            if sample in treated_samples_full[breaking]['Pan']['YES']:
                forbidden = True
        if forbidden is not True:
            not_breaked.add(sample)

    return only_treated, not_breaked


def do_plot(samples_lowbraca, samples_braca, not_breaked, only_treated, exp):

    dic_primary_full, _ = return_metadata()
    color_ttype = return_colors()
    fig, ax = plt.subplots(1, 1, figsize=(2, 1.5))
    config_params(6.5)

    sigs = ['12_ID6_0.962941_1']

    exposures_not_breaked_1 = exp[sigs].sum(axis=1).loc[[i for i in not_breaked if i in samples_lowbraca]].dropna()
    exposures_not_breaked_2 = exp[sigs].sum(axis=1).loc[[i for i in not_breaked if i in samples_braca]].dropna()
    exposures_breaked_1 = exp[sigs].sum(axis=1).loc[[i for i in only_treated if i in samples_lowbraca]].dropna()
    exposures_breaked_2 = exp[sigs].sum(axis=1).loc[[i for i in only_treated if i in samples_braca]].dropna()

    sns.boxplot(
        data=[
            exposures_not_breaked_1, exposures_breaked_1, exposures_not_breaked_2, exposures_breaked_2
        ],
        linewidth=0.6,
        showfliers=False,
        color='#cbcacbff'
    )

    plt.ylabel('Indels DSB repair by\nnon-homologous end-joining')
    plt.xticks(
        [0, 1, 2, 3],
        [
            'Not radiated no BRCAness ({})'.format(len(exposures_not_breaked_1)),
            'Radiated no BRCAness ({})'.format(len(exposures_breaked_1)),
            'Not radiated BRCAness ({})'.format(len(exposures_not_breaked_2)),
            'Radiated BRCAness ({})'.format(len(exposures_breaked_2))
        ],
        rotation=90
    )

    plotdot = []
    colors = []
    for sample in [i for i in not_breaked if i in samples_lowbraca]:
        plotdot.append(exp[sigs].sum(axis=1).loc[sample])
        colors.append(color_ttype[dic_primary_full[sample]])

    ax.scatter(
        [0 + np.random.uniform(-0.2, 0.2, 1)[0] for i in range(len(plotdot))],
        plotdot,
        color=colors,
        s=1,
        alpha=0.2
    )

    plotdot = []
    colors = []
    for sample in [i for i in only_treated if i in samples_lowbraca]:
        plotdot.append(exp[sigs].sum(axis=1).loc[sample])
        colors.append(color_ttype[dic_primary_full[sample]])

    ax.scatter(
        [1 + np.random.uniform(-0.2, 0.2, 1)[0] for i in range(len(plotdot))],
        plotdot,
        color=colors,
        s=1,
        alpha=0.2
    )

    plotdot = []
    colors = []
    for sample in [i for i in not_breaked if i in samples_braca]:
        plotdot.append(exp[sigs].sum(axis=1).loc[sample])
        colors.append(color_ttype[dic_primary_full[sample]])

    ax.scatter(
        [2 + np.random.uniform(-0.2, 0.2, 1)[0] for i in range(len(plotdot))],
        plotdot,
        color=colors,
        s=1,
        alpha=0.2
    )

    plotdot = []
    colors = []
    for sample in [i for i in only_treated if i in samples_braca]:
        plotdot.append(exp[sigs].sum(axis=1).loc[sample])
        colors.append(color_ttype[dic_primary_full[sample]])

    ax.scatter(
        [3 + np.random.uniform(-0.2, 0.2, 1)[0] for i in range(len(plotdot))],
        plotdot,
        color=colors,
        s=1,
        alpha=0.2
    )

    plt.ylim(0, 700)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig('figures/radiation.svg')
    plt.show()

    ##################

    fig, ax = plt.subplots(1, 1, )
    stat, pval1 = mannwhitneyu(exposures_not_breaked_1, exposures_breaked_1)
    print("Not radiated no BRCAnes vs Radiated no BRCAness", pval1)

    stat, pval2 = mannwhitneyu(exposures_not_breaked_2, exposures_breaked_2)
    print("Not radiated  BRCAnes vs Radiated  BRCAness", pval2)

    stat, pval3 = mannwhitneyu(exposures_breaked_1, exposures_breaked_2)
    print("radiated  no BRCAnes vs Radiated  BRCAness", pval3)

    ax.text(1, 1, "$\it{P}$" + " = {}".format(sci_notation(pval1)), fontsize=7)
    ax.text(1, 4, "$\it{P}$" + " = {}".format(sci_notation(pval2)), fontsize=7)
    ax.text(1, 8, "$\it{P}$" + " = {}".format(sci_notation(pval3)), fontsize=7)

    plt.xlim(0, 5)
    plt.ylim(0, 10)
    plt.savefig('figures/radiation_pvals.svg')
    sys.exit()


def run():

    config_params()
    np.random.seed(12345)
    treated_samples_full = pickle.load(gzip.open(Conf['treatment_FDA_drug']))

    base_path = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/'
    exposures_path = base_path + 'snvs/exposures/Pan_full/Pan_full.exposures.tsv'
    # load exposures of Pan
    exp_snv = load_exposures(exposures_path)

    exposures_path = base_path + 'indels/exposures/Pan/Pan.exposures.tsv'

    # load exposures of Pan
    exp = load_exposures(exposures_path)

    treatment = 'RADIATION'

    breaks = [
        'Unknown', 'Topoisomerase Inhibitor', 'Anthracycline Topoisomerase Inhibitor',
        'Alkylating Drug', 'Platinum-based Drug', 'Poly(ADP-Ribose) Polymerase Inhibitor',
        'Miscellanious', 'Nucleoside Analog Antiviral', 'TOPOII', 'Nuclear therapy'
    ]

    samples_braca, samples_lowbraca = select_HRdeficient_samples(exp_snv)
    only_treated, not_breaked = select_treated_not_treated(treated_samples_full, treatment, breaks, exp)

    do_plot(samples_lowbraca, samples_braca, not_breaked, only_treated, exp)


if __name__ == '__main__':
    run()
