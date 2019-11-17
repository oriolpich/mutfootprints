# Import modules
import sys
import os
from collections import defaultdict
import pickle
import gzip

import pandas as pd
from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from utils.metadata import return_metadata
from utils.colors_ttype import return_colors
from utils.plots_utils import config_params


# form stackoverflow
def sci_notation(number, sig_fig=2):

    ret_string = "{0:.{1:d}e}".format(number, sig_fig)
    a, b = ret_string.split("e")
    b = int(b)  # removed leading "+" and strips leading zeros too.
    string = a + r"x$10^{" + str(b) + "}$"

    return string


def get_distribution_treatment(drug, file_days_treatment, exposures_path):

    counter_length = pickle.load(gzip.open(file_days_treatment))
    dic_ttypes, _ = return_metadata()
    _, samples_exp = read_exposures(exposures_path)

    keep_samples_d = defaultdict(list)
    days_treat_dict = defaultdict(list)

    for sample, v in counter_length.items():
        if drug in v:
            numbers = [i for i in v[drug] if str(i) != 'nan']
            if len(numbers) > 0:
                if str(dic_ttypes[sample]) != 'nan':
                    days = np.sum(numbers)
                    days_treat_dict[dic_ttypes[sample]].append(days)
                    keep_samples_d[dic_ttypes[sample]].append(sample)

    # get the mean of each tumor type
    keep_mean = {}
    for ttype, l in days_treat_dict.items():
        if (str(ttype) != 'nan') & (str(ttype) != 'Unknown'):
            if len(l) >= 10:
                keep_mean[ttype] = np.mean(l)

    ttype_ordereds = [i for i in list(sorted(keep_mean, key=keep_mean.get))]

    return ttype_ordereds, days_treat_dict, keep_samples_d


def plot_long_short_exposure(drug, exposures_path,  signatures_snv, ttype_ordereds, days_treat_dict, keep_samples_d, ymax, type_mut, type_method, sig1 = False):

    dchange = {'snv': 'SBS', 'dbs': 'DBS'}
    fix_type_mut = dchange[type_mut]
    all_sigs = signatures_snv[drug]
    df_exposures, _ = read_exposures(exposures_path)
    color_ttype = return_colors()

    tobox = []
    colors = []
    samples_exp = df_exposures.index.tolist()
    labels = []
    keep_pvals = []

    for ttype in ttype_ordereds:

        percentile_long = np.percentile(days_treat_dict[ttype], 75)
        sample_long = [keep_samples_d[ttype][ix] for ix, s in enumerate(days_treat_dict[ttype]) if s >= percentile_long]
        percentile_low = np.percentile(days_treat_dict[ttype], 25)

        sample_low = [keep_samples_d[ttype][ix] for ix, s in enumerate(days_treat_dict[ttype]) if s <= percentile_low]
        samplelongind = [s for s in sample_long if s in samples_exp]
        samplelowind = [s for s in sample_low if s in samples_exp]

        total_number_samples = [keep_samples_d[ttype][ix] for ix, s in enumerate(days_treat_dict[ttype])]
        if (len(samplelongind) > 25) & (len(samplelowind) > 25):

            if (df_exposures.loc[samplelongind][all_sigs].sum().sum() > 0) & (df_exposures.loc[samplelowind][all_sigs].sum().sum() > 0):
                list_short = df_exposures.loc[samplelowind][all_sigs].sum(axis=1)
                list_long = df_exposures.loc[samplelongind][all_sigs].sum(axis=1)

                tobox.append(list_short)
                tobox.append(list_long)
                stat, pval = mannwhitneyu(list_short, list_long)
                keep_pvals.append(pval)
                colors.append(color_ttype[ttype])
                colors.append(color_ttype[ttype])
                labels.append('{} ST\n({})'.format(ttype, len(list_short)))
                labels.append('{} LT\n({})'.format(ttype, len(list_long)))

    if len(tobox) == 2:

        config_params(5.5)
        fig, ax = plt.subplots(1, 1, figsize=(len(tobox) / 3, 1))
        box = sns.boxplot(data=tobox,  color='black', boxprops=dict(alpha=.9, ), linewidth=0.6, showfliers=False,) #color='#ccccccff')

        for b in box.artists:
            b.set_facecolor('#d1d1d1ff')

        ax.scatter([0 + np.random.uniform(-0.2, 0.2, 1)[0] for i in tobox[0]],
                   tobox[0], color=colors[0], s=1, alpha=0.6)
        ax.scatter([1 + np.random.uniform(-0.2, 0.2, 1)[0] for i in tobox[1]], tobox[1], color=colors[1], s=1,
                   alpha=0.6)
        ax.set_xticklabels(labels, rotation=90)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.title("$\it{P}$" + " = {}".format(sci_notation(pval)), fontsize=5)

        plt.ylim(0, ymax)

        if sig1 is False:
            plt.ylabel('{} related {}'.format(drug.lower().capitalize(), fix_type_mut))
        else:
            plt.ylabel('Aging-related {}'.format(fix_type_mut))

        if sig1 is False:
            plt.savefig(
                'figures/{}/length_treatment_exposure_{}_{}.svg'.format(type_method, drug, type_mut)
            )
            plt.savefig(
                'figures/{}/length_treatment_exposure_{}_{}.png'.format(type_method, drug, type_mut),
                dpi=600, bbox_inches="tight"
            )
        else:
            plt.savefig(
                'figures/{}/length_treatment_exposure_{}_{}_sig1.svg'.format(type_method, drug, type_mut)
            )
            plt.savefig(
                'figures/{}/length_treatment_exposure_{}_{}_sig1.png'.format(type_method, drug, type_mut),
                dpi=600, bbox_inches="tight"
            )
        plt.close()

    else:
        fig, axs = plt.subplots(1, 2, figsize=(1.25, 1))

        axs[0].set_title("$\it{P}$" + " = {}".format(sci_notation(keep_pvals[0])), fontsize=5)
        axs[1].set_title("$\it{P}$" + " = {}".format(sci_notation(keep_pvals[1])), fontsize=5)

        box = sns.boxplot(
            data=tobox[0:2],  color='black', boxprops=dict(alpha=.9, ),
            linewidth=0.6, showfliers=False, ax=axs[0]
        )  # color='#ccccccff')

        for b in box.artists:
            b.set_facecolor('#d1d1d1ff')

        box = sns.boxplot(
            data=tobox[2:], color='black', boxprops=dict(alpha=.9, ),
            linewidth=0.6, showfliers=False, ax=axs[1]
        )  # color='#ccccccff')

        for b in box.artists:
            b.set_facecolor('#d1d1d1ff')

        axs[0].scatter(
            [0 + np.random.uniform(-0.2, 0.2, 1)[0] for i in tobox[0]],
            tobox[0], color=colors[0], s=1, alpha=0.6
        )
        axs[0].scatter(
            [1 + np.random.uniform(-0.2, 0.2, 1)[0] for i in tobox[1]],
            tobox[1], color=colors[1], s=1, alpha=0.6
        )
        axs[1].scatter(
            [0 + np.random.uniform(-0.2, 0.2, 1)[0] for i in tobox[2]],
            tobox[2], color=colors[2], s=1, alpha=0.6
        )
        axs[1].scatter(
            [1 + np.random.uniform(-0.2, 0.2, 1)[0] for i in tobox[3]],
            tobox[3], color=colors[3], s=1, alpha=0.6
        )
        axs[0].spines['top'].set_visible(False)
        axs[0].spines['right'].set_visible(False)
        axs[1].spines['top'].set_visible(False)
        axs[1].spines['right'].set_visible(False)

        axs[0].set_xticklabels(labels[:2], rotation=90)
        axs[1].set_xticklabels(labels[2:], rotation=90)

        axs[0].set_ylabel('{} related {}'.format(drug.lower().capitalize(), fix_type_mut))
        axs[1].set_ylabel('{} related {}'.format(drug.lower().capitalize(), fix_type_mut))
        axs[0].set_ylim(0, ymax)

        axs[1].set_ylim(0, ymax)

        if sig1 is False:
            plt.ylabel('{} related {}'.format(drug.lower().capitalize(), fix_type_mut))
        else:
            plt.ylabel('Aging-related {}'.format(fix_type_mut))

        if sig1 is False:
            plt.savefig(
                'figures/{}/length_treatment_exposure_{}_{}.svg'.format(type_method, drug, type_mut)
            )
            plt.savefig(
                'figures/{}/length_treatment_exposure_{}_{}.png'.format(type_method, drug, type_mut),
                dpi=600, bbox_inches="tight"
            )
        else:
            plt.savefig(
                'figures/{}/length_treatment_exposure_{}_{}_sig1.svg'.format(type_method, drug, type_mut)
            )
            plt.savefig(
                'figures/{}/length_treatment_exposure_{}_{}_sig1.png'.format(type_method, drug, type_mut),
                dpi=600, bbox_inches="tight"
            )
        plt.close()


def read_exposures(exposures_path):

    df = pd.read_csv(exposures_path, sep='\t', index_col=0)
    df = df.T
    samples_exp = df.index.tolist()
    return df, samples_exp


def run():

    np.random.seed(1234)
    # ---------------
    # Necessary files
    # --------------
    # file with days of treatment calculations
    file_days_treatment = 'data/clinical_data/dic_counter_total_days_treatment.pckl.gz'

    # ------------
    # SigProfiler
    # ------------
    signatures_snv = {
        'CISPLATIN': ['5_DBS5_0.944431_1.0'],
        'CARBOPLATIN': ['5_DBS5_0.944431_1.0'],
        'OXALIPLATIN': ['5_DBS5_0.944431_1.0'],
    }
    type_method = 'SigProfiler'

    print('SigProfiler DBS...')

    # file with exposures
    exposures_path = 'data/hartwig/signatures/extraction/results/SigProfiler/dbs/exposures/'
    exposures_path += 'PanNoSkinNoUnknown.dbs/PanNoSkinNoUnknown.dbs.exposures.tsv'

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CISPLATIN', file_days_treatment, exposures_path
    )
    ylim = 500
    plot_long_short_exposure(
        'CISPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "dbs", type_method
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CARBOPLATIN', file_days_treatment, exposures_path
    )
    ylim = 300
    plot_long_short_exposure(
        'CARBOPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "dbs", type_method
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'OXALIPLATIN', file_days_treatment, exposures_path
    )
    ylim = 500
    plot_long_short_exposure(
        'OXALIPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "dbs", type_method
    )

    signatures_snv = {
        'CISPLATIN': ['1_SBS31_0.968153_0.98'],
        'CARBOPLATIN': ['1_SBS31_0.968153_0.98'],
        'CAPECITABINE': ['19_SBS17b_0.961548_0.99'],
        '5-FU_CAPE': ['19_SBS17b_0.961548_0.99'],
        'OXALIPLATIN': ['20_0.92']
    }

    # file with exposures
    exposures_path = 'data/hartwig/signatures/extraction/results/SigProfiler/snvs/exposures/'
    exposures_path += 'PanNoSkinNoUnknown.snvs/PanNoSkinNoUnknown.snvs.exposures.tsv'

    print('SigProfiler SNVs...')
    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        '5-FU_CAPE', file_days_treatment, exposures_path
    )

    print(ttype_ordereds, keep_samples_d)
    ylim = 11000
    plot_long_short_exposure(
        '5-FU_CAPE',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CISPLATIN', file_days_treatment, exposures_path
    )
    ylim = 11000
    plot_long_short_exposure(
        'CISPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CARBOPLATIN', file_days_treatment, exposures_path
    )
    ylim = 4000
    plot_long_short_exposure(
        'CARBOPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'OXALIPLATIN', file_days_treatment, exposures_path
    )
    ylim = 10000
    plot_long_short_exposure(
        'OXALIPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CAPECITABINE', file_days_treatment, exposures_path
    )
    ylim = 10000
    plot_long_short_exposure(
        'CAPECITABINE',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method
    )

    # SIGNATURE 1

    signatures_snv = {
        'CISPLATIN': ['23_SBS1_0.996494_1.0'],
        'CARBOPLATIN': ['23_SBS1_0.996494_1.0'],
        '5-FU_CAPE': ['23_SBS1_0.996494_1.0'],
        'OXALIPLATIN': ['23_SBS1_0.996494_1.0']
    }

    print('SigProfiler SBS1 SNVs...')

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CISPLATIN', file_days_treatment, exposures_path
    )
    ylim = 11000
    plot_long_short_exposure(
        'CISPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method, sig1=True
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CARBOPLATIN', file_days_treatment, exposures_path
    )
    ylim = 4000
    plot_long_short_exposure(
        'CARBOPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method, sig1=True
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'OXALIPLATIN', file_days_treatment, exposures_path
    )
    ylim = 10000

    plot_long_short_exposure(
        'OXALIPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method, sig1=True
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        '5-FU_CAPE', file_days_treatment, exposures_path
    )
    ylim = 10000
    plot_long_short_exposure(
        '5-FU_CAPE',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method, sig1=True
    )

    # =============
    # SignatureAnalyzer
    # =============
    
    print('SignatureAnalyzer DBS...')

    signatures_dbs = {
        'CISPLATIN': ['9_1', '3_1'],
        'CARBOPLATIN': ['9_1', '3_1'],
        'OXALIPLATIN': ['9_1', '3_1']
    }

    # file with exposures
    exposures_path = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/'
    exposures_path += 'dbs/exposures/Pan/Pan.exposures.tsv'

    type_method = 'SignatureAnalyzer'
    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CISPLATIN', file_days_treatment, exposures_path
    )
    ylim = 500
    plot_long_short_exposure(
        'CISPLATIN', exposures_path,  signatures_dbs, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "dbs", type_method
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CARBOPLATIN', file_days_treatment, exposures_path
    )
    ylim = 200
    plot_long_short_exposure(
        'CARBOPLATIN', exposures_path,  signatures_dbs, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "dbs", type_method
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'OXALIPLATIN', file_days_treatment, exposures_path
    )
    ylim = 500
    plot_long_short_exposure(
        'OXALIPLATIN', exposures_path,  signatures_dbs, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "dbs", type_method
    )

    print('SignatureAnalyzer SNVs...')

    # file with exposures
    exposures_path = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/'
    exposures_path += 'snvs/exposures/Pan_full/Pan_full.exposures.tsv'

    signatures_snv = {
        'CISPLATIN': ['21_SBS31_0.953955_1', '14_1'],
        'CARBOPLATIN': ['21_SBS31_0.953955_1', '25_1'],
        '5-FU_CAPE': ['31_SBS17b_0.968799_1'],
        'OXALIPLATIN': ['14_1', '37_1']
    }
    # -------------
    # Run for each treatment
    # -------------

    type_method = 'SignatureAnalyzer'
    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CISPLATIN', file_days_treatment, exposures_path
    )
    ylim = 6000
    plot_long_short_exposure(
        'CISPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CARBOPLATIN', file_days_treatment, exposures_path
    )
    ylim = 2000
    plot_long_short_exposure(
        'CARBOPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'OXALIPLATIN', file_days_treatment, exposures_path
    )
    ylim = 6000
    plot_long_short_exposure(
        'OXALIPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        '5-FU_CAPE', file_days_treatment, exposures_path
    )
    ylim = 10000
    plot_long_short_exposure(
        '5-FU_CAPE',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method
    )

    signatures_snv = {
        'CISPLATIN': ['17_SBS1_0.995019_1'],
        'CARBOPLATIN': ['17_SBS1_0.995019_1'],
        '5-FU_CAPE': ['17_SBS1_0.995019_1'],
        'OXALIPLATIN': ['17_SBS1_0.995019_1']
    }

    print('SignatureAnalyzer SBS1 SNVs...')

    type_method = 'SignatureAnalyzer'
    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CISPLATIN', file_days_treatment, exposures_path
    )
    ylim = 6000
    plot_long_short_exposure(
        'CISPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method, sig1=True
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'CARBOPLATIN', file_days_treatment, exposures_path
    )
    ylim = 2000
    plot_long_short_exposure(
        'CARBOPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method, sig1=True
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        'OXALIPLATIN', file_days_treatment, exposures_path
    )
    ylim = 6000
    plot_long_short_exposure(
        'OXALIPLATIN',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method, sig1=True
    )

    ttype_ordereds, days_treat_dict, keep_samples_d = get_distribution_treatment(
        '5-FU_CAPE', file_days_treatment, exposures_path
    )
    ylim = 6000
    plot_long_short_exposure(
        '5-FU_CAPE',  exposures_path, signatures_snv, ttype_ordereds,
        days_treat_dict, keep_samples_d, ylim, "snv", type_method, sig1=True
    )


if __name__ == '__main__':
    run()
