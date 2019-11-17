# Import module
import os
import sys
import gzip
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import dill as pickle
from tqdm import tqdm

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..',))
sys.path.insert(0, scripts_path)

from config import Conf
from utils.metadata import return_metadata
from utils.colors_ttype import return_colors
from utils.plots_utils import config_params


def calculate_proportions_timing(file_snvs, file_dbs, file_indels, signature_2_treatment):

    treated_samples_dic = pickle.load(gzip.open(Conf['treatment_specific_drug']))
    outfile_all_cons = '{}/clonality_muts.gz'.format(os.path.dirname(file_snvs))

    # read files and get only selected signatures
    if not os.path.isfile(outfile_all_cons):
        print('loading SNV')

        df_snv = pd.read_csv(file_snvs, sep='\t', usecols=['SAMPLE', 'ML', 'timing_class'], low_memory=False)
        df_snv['ML'] = df_snv['ML'].apply(lambda x: '{}_snvs'.format(x))
        df_snv = df_snv[df_snv['ML'].isin(signature_2_treatment)]

        print('loading DBS')
        df_dbs = pd.read_csv(file_dbs, sep='\t', usecols=['SAMPLE', 'ML', 'timing_class'], low_memory=False)
        df_dbs['ML'] = df_dbs['ML'].apply(lambda x: '{}_dbs'.format(x))
        df_dbs = df_dbs[df_dbs['ML'].isin(signature_2_treatment)]

        print('loading ID')
        df_indels = pd.read_csv(file_indels, sep='\t', usecols=['SAMPLE', 'ML', 'timing_class'], low_memory=False)
        df_indels['ML'] = df_indels['ML'].apply(lambda x: '{}_indels'.format(x))
        df_indels = df_indels[df_indels['ML'].isin(signature_2_treatment)]

        total_merged = pd.concat([df_snv, df_dbs, df_indels])
        total_merged.to_csv(outfile_all_cons, sep='\t', index=False, header=True)
    else:
        total_merged = pd.read_csv(outfile_all_cons, sep='\t')

    d_box = defaultdict(list)
    dttype = defaultdict(lambda: defaultdict(list))
    d_box_subclonality = defaultdict(list)
    dttype_subclonality = defaultdict(lambda: defaultdict(list))
    dic_primary_full, _ = return_metadata()

    for sample, data in tqdm(total_merged.groupby(by='SAMPLE')):

        for sig in signature_2_treatment:

            affected = False

            # select only samples treated with the drug
            for drug in signature_2_treatment[sig]:
                set_samples = treated_samples_dic[drug]['Pan']['YES']
                if sample in set_samples:
                    affected = True

            if affected is True:
                # amount of mutations of each type per sample
                total_early = len(data[data['timing_class'] == 'clonal [early]'])
                total_late = len(data[data['timing_class'] == 'clonal [late]'])
                total_early_sum = total_early + total_late
                total_subclonal = len(data[data['timing_class'] == 'subclonal'])

                # select the cases of the signature we are interested in
                filtered = data[data['ML'] == sig]

                # calculate the different activities
                if len(filtered):
                    nclonal = len(filtered[filtered['timing_class'] == 'clonal [early]'])
                    nlate = len(filtered[filtered['timing_class'] == 'clonal [late]'])
                    n_early_sum = nclonal + nlate
                    nsubclonal = len(filtered[filtered['timing_class'] == 'subclonal'])

                    # proportion clonal early / clonal late
                    if (nclonal > 5) & (nlate > 5):
                        proportion_late = nlate / total_late
                        proportion_early = nclonal / total_early
                        d_box[sig].append(proportion_late / proportion_early)
                        dttype[sig][dic_primary_full[sample]].append(proportion_late / proportion_early)

                    # proportion subclonal/other
                    if (n_early_sum > 5) & (nsubclonal > 5):
                        proportion_subclonal = nsubclonal / total_subclonal
                        proportion_early = n_early_sum / total_early_sum
                        d_box_subclonality[sig].append(proportion_subclonal / proportion_early)
                        dttype_subclonality[sig][dic_primary_full[sample]].append(proportion_subclonal / proportion_early)

    return d_box, dttype, d_box_subclonality, dttype_subclonality


def get_values(file_snvs, file_dbs, file_indels, signature_2_treatment):

    outfile_clon = '{}/clonality_dict.pckl.gz'.format(os.path.isfile(file_snvs))
    outfile_subclon = '{}/subclonality_dict.pckl.gz'.format(os.path.isfile(file_snvs))

    if (os.path.isfile(outfile_clon)) & (os.path.isfile(outfile_subclon)):
        d_box = pickle.load(gzip.open(outfile_clon))
        dttype = pickle.load(gzip.open(outfile_clon.replace('dict', 'dict_ttype')))
        d_box_subclonality = pickle.load(gzip.open(outfile_clon))
        dttype_subclonality = pickle.load(gzip.open(outfile_subclon.replace('dict', 'dict_ttype')))

    else:
        d_box, dttype, d_box_subclonality, dttype_subclonality =  calculate_proportions_timing(file_snvs, file_dbs, file_indels, signature_2_treatment)

    return d_box, dttype, d_box_subclonality, dttype_subclonality


def do_plot(d_box, dttype, outfile, signature_2_treatment):

    np.random.seed(12345)

    config_params(6)

    toplot = []
    labels = []

    for sig in sorted(d_box, key=lambda k: np.median(d_box[k]), reverse=False):

        if (len(d_box[sig]) > 5) & (sig in signature_2_treatment):
            toplot.append(sorted(d_box[sig]))
            labels.append(sig)

    fig, ax = plt.subplots(1, 1, figsize=(3, 3))

    plt.yscale('log')
    color_ttype = return_colors()

    number_of_samples = []
    for ix, signature in enumerate(labels):
        plotdot = []
        colors = []
        for ttype, samples in dttype[signature].items():
            for sample in samples:
                plotdot.append(sample)
                colors.append(color_ttype[ttype])
        ax.scatter([ix + np.random.uniform(-0.2, 0.2, 1)[0] for i in range(len(plotdot))],
                   plotdot, color=colors, s=1, alpha=0.75)
        number_of_samples.append(len(plotdot))

    box = sns.boxplot(data=toplot, color='black', linewidth=0.6, ax=ax, showfliers=False)
    for b in box.artists:
        b.set_facecolor('#d1d1d1ff')

    plt.ylabel('foldchange late/early')
    ax.set_xticklabels(['{}\n{}'.format(l, number_of_samples[ixl]) for ixl, l in enumerate(labels)], rotation=90)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0.05, 25)

    ax.hlines(1, 0, len(toplot), alpha=0.5)

    plt.tight_layout()
    plt.savefig(outfile)


def main(file_snvs, file_dbs, file_indels, signature_2_treatment, type_extraction):

    d_box, dttype, d_box_subclonality, dttype_subclonality = get_values(
        file_snvs, file_dbs, file_indels, signature_2_treatment
    )

    os.makedirs('figures/{}'.format(type_extraction), exist_ok=True)
    do_plot(d_box, dttype, 'figures/{}/early_late_fig.svg'.format(type_extraction), signature_2_treatment)
    do_plot(d_box_subclonality, dttype_subclonality, 'figures/{}/subclonal_clonal_fig.svg'.format(type_extraction), signature_2_treatment)


def plot_timing_analysis():

    file_snvs = 'data/hartwig/timing/merged/ML/SignatureAnalyzer/Pan.ML.snvs.gz'
    file_dbs = 'data/hartwig/timing/merged/ML/SignatureAnalyzer/Pan.ML.dbs.gz'
    file_indels = 'data/hartwig/timing/merged/ML/SignatureAnalyzer/Pan.ML.indels.gz'

    type_extraction = 'SignatureAnalyzer'

    signature_2_treatment = {
        '9_1_dbs': ['OXALIPLATIN', 'CISPLATIN', 'CARBOPLATIN'],
        '3_1_dbs': ['OXALIPLATIN', 'CISPLATIN', 'CARBOPLATIN'],
        '37_1_snvs': ['OXALIPLATIN'],
        '14_1_snvs': ['OXALIPLATIN', 'CISPLATIN'],

        '21_SBS31_0.953955_1_snvs': ['CISPLATIN', 'CARBOPLATIN'],
        '31_SBS17b_0.968799_1_snvs': ['5-FU_CAPE'],
        '25_1_snvs': ['CARBOPLATIN'],

        '17_SBS1_0.995019_1_snvs': ['CISPLATIN', 'CARBOPLATIN', 'OXALIPLATIN', '5-FU_CAPE'],
        '6_SBS4_0.961549_1_snvs': ['CISPLATIN', 'CARBOPLATIN', 'OXALIPLATIN', '5-FU_CAPE']
    }

    main(file_snvs, file_dbs, file_indels, signature_2_treatment, type_extraction)

    # ----------------
    # SigProfiler
    # ----------------
    file_snvs = 'data/hartwig/timing/merged/ML/SigProfiler/Pan.ML.snvs.gz'
    file_dbs = 'data/hartwig/timing/merged/ML/SigProfiler/Pan.ML.dbs.gz'
    file_indels = 'data/hartwig/timing/merged/ML/SigProfiler/Pan.ML.indels.gz'
    type_extraction = 'SigProfiler'

    signature_2_treatment = {
        '20_0.92_snvs': ['OXALIPLATIN'],
        '5_DBS5_0.944431_1.0_dbs': ['CARBOPLATIN', 'CISPLATIN', 'OXALIPLATIN'],
        '19_SBS17b_0.961548_0.99_snvs': ['5-FU_CAPE'],
        '1_SBS31_0.968153_0.98_snvs': ['CISPLATIN', 'CARBOPLATIN'],
        '23_SBS1_0.996494_1.0_snvs': ['CISPLATIN', 'CARBOPLATIN', 'OXALIPLATIN', 'CAPECITABINE'],
        '17_SBS4_0.935484_0.97_snvs': ['CISPLATIN', 'CARBOPLATIN', 'OXALIPLATIN', 'CAPECITABINE'],
    }

    main(file_snvs, file_dbs, file_indels, signature_2_treatment, type_extraction)


if __name__ == '__main__':
    plot_timing_analysis()
