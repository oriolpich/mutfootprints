# Import modules
import sys
import os
from collections import defaultdict
import gzip

import pandas as pd
import dill as pickle
from scipy.stats import pearsonr
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from utils.metadata import return_metadata
from utils.plots_utils import config_params
from config import Conf


# Fix the random seed
np.random.seed(12345)


def read_exposures(exposures_path):
    df = pd.read_csv(exposures_path, sep='\t', index_col=0)
    df = df.T
    return df


def plot_ttype_distribution(dsigprof, dsiganalyzer):

    for ttype in dsigprof:
        samples_sorted = sorted(dsigprof[ttype], key=dsigprof[ttype].get)
        xvals = []
        yvals = []
        tots = []
        fig, ax = plt.subplots(1, 1, figsize=(2, 2.5))

        for ix, s in enumerate(samples_sorted):
            if (dsigprof[ttype][s] > 0) & (dsiganalyzer[ttype][s] > 0):
                vals = [dsigprof[ttype][s], dsiganalyzer[ttype][s]]
                plt.vlines(ix, np.min(vals), np.max(vals))
                plt.scatter(ix, vals[0], c='darkred', s=2)
                plt.scatter(ix, vals[1], c='darkblue', s=2)

                yvals.append(vals[0])
                xvals.append(vals[1])
                tots.append(abs(vals[0] - vals[1]))

        plt.title(ttype)
        plt.close()


def plot_single_distribution(dsigprof, dsiganalyzer, drug, outpath, type_dic):
    samples_sorted = sorted(dsigprof, key=dsigprof.get)
    xvals = []
    yvals = []
    fig, ax = plt.subplots(1, 1, figsize=(2, 1))

    for ix, s in enumerate(samples_sorted):
        if (dsigprof[s] > 0) & (dsiganalyzer[s] > 0):
            vals = [dsigprof[s], dsiganalyzer[s]]
            plt.vlines(ix, np.min(vals), np.max(vals), lw=0.4, color='grey', alpha=0.4)
            plt.scatter(ix, vals[0], c='darkred', s=2.5)
            plt.scatter(ix, vals[1], c='darkblue', s=2.5)

            yvals.append(vals[0])
            xvals.append(vals[1])

    stat, pval = pearsonr(xvals, yvals)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Exposure')
    ax.set_xlabel('Samples')
    plt.title('{} pearson = {}, pval = {}, n = {}'.format(drug, round(stat, 3), round(pval, 6), len(xvals)))
    plt.savefig('{}/{}_{}.svg'.format(outpath, drug, type_dic))
    plt.close()


def get_samples_similar_exposure(dsigprof, dsiganalyzer, drug, outpath, type_dic):
    samples_sorted = sorted(dsigprof, key=dsigprof.get)
    fig, ax = plt.subplots(1, 1, figsize=(1.5, 1.5))
    xvals = []
    yvals = []
    bad_samples = set()

    for ix, s in enumerate(samples_sorted):

        # check if we have exposure in both cases, otherwise they should go to a different category
        if (dsigprof[s] > 0) & (dsiganalyzer[s] > 0):

            vals = [dsigprof[s], dsiganalyzer[s]]
            distance = abs(vals[0] - vals[1])
            if distance <= 0.15:

                plt.vlines(ix, np.min(vals), np.max(vals))
                plt.scatter(ix, vals[0], c='darkred', s=2)
                plt.scatter(ix, vals[1], c='darkblue', s=2)

                yvals.append(vals[0])
                xvals.append(vals[1])
            else:
                bad_samples.add(s)
        else:
            bad_samples.add(s)
    stat, pval = pearsonr(xvals, yvals)
    plt.title('{} pearson = {}, pval = {}, n = {}'.format(drug, round(stat, 3), round(pval, 6), len(xvals)))
    plt.savefig('{}/{}_{}.coherent.png'.format(outpath, drug, type_dic), dpi=600, bbox_inches='tight')
    plt.close()

    return bad_samples


def barplot_classes(dsigprof, dsiganalyzer, drug, outpath, type_dic):
    samples_sorted = sorted(dsigprof, key=dsigprof.get)

    fig, ax = plt.subplots(1, 1, figsize=(1, 1))

    class_sigpro = 0
    class_sigan = 0
    class_0 = 0
    class_disc = 0
    class_conc = 0

    for ix, s in enumerate(samples_sorted):

        if (dsigprof[s] == 0) & (dsiganalyzer[s] < 0.0075):
            class_0 += 1
        elif (dsigprof[s] > 0) & (dsiganalyzer[s] < 0.0075):
            class_sigpro += 1
        elif (dsigprof[s] == 0) & (dsiganalyzer[s] >= 0.0075):
            class_sigan += 1
        elif (dsigprof[s] > 0) & (dsiganalyzer[s] >= 0.0075):

            vals = [dsigprof[s], dsiganalyzer[s]]
            distance = abs(vals[0] - vals[1])

            if distance < 0.15:
                class_conc += 1
            else:
                class_disc += 1
    bottom = 0

    for v in [class_conc, class_disc, class_sigan, class_sigpro, class_0]:
        val = v / len(samples_sorted)
        plt.bar(0, val, bottom=bottom)
        bottom += val
    plt.savefig('{}/{}_{}.barplot.png'.format(outpath, drug, type_dic), dpi=600, bbox_inches='tight')
    plt.close()


def plot_single_correlation(dsigprof, dsiganalyzer, drug, outpath, type_dic):
    samples_sorted = sorted(dsigprof, key=dsigprof.get)
    xvals = []
    yvals = []
    fig, ax = plt.subplots(1, 1, figsize=(1, 1))

    for ix, s in enumerate(samples_sorted):
        if (dsigprof[s] > 0) & (dsiganalyzer[s] > 0):
            vals = [dsigprof[s], dsiganalyzer[s]]
            plt.scatter(dsigprof[s], dsiganalyzer[s], c='grey', s=1)
            yvals.append(vals[0])
            xvals.append(vals[1])

    stat, pval = pearsonr(xvals, yvals)
    plt.ylabel('Signature Analyzer')
    plt.xlabel('SigProfiler')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.title('{} pearson = {}, pval = {}, n = {}'.format(drug, round(stat, 3), round(pval, 6), len(xvals)))
    plt.savefig('{}/{}_{}.correlation.svg'.format(outpath, drug, type_dic))

    if type_dic == 'exposure':

        fig, ax = plt.subplots(1, 1, figsize=(1, 1))
        xvals = []
        yvals = []

        for ix, s in enumerate(samples_sorted):
            if (dsigprof[s] > 0) & (dsiganalyzer[s] > 0):
                vals = [dsigprof[s], dsiganalyzer[s]]

                distance = abs(vals[0] - vals[1])
                if distance <= 0.15:
                    plt.scatter(dsigprof[s], dsiganalyzer[s], c='grey', s=1)
                    yvals.append(vals[0])
                    xvals.append(vals[1])

        stat, pval = pearsonr(xvals, yvals)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel('SigProfiler')
        ax.set_ylabel('Signature Analyzer')

        plt.title('{} pearson = {}, pval = {}, n = {}'.format(drug, round(stat, 3), round(pval, 6), len(xvals)))
        plt.savefig('{}/{}_{}.correlation.coherent.svg'.format(outpath, drug, type_dic))
        plt.close()


def process_exposures_methods(
        path_siganalyzer, path_sigprofiler, drug,
        signature_2_treatment_signature_analyzer,
        signature_2_treatment_sigprofiler, outpath):

    df_sigpro = read_exposures(path_sigprofiler)
    df_sigan = read_exposures(path_siganalyzer)
    all_treated_samples = set()

    treated_samples = pickle.load(gzip.open(Conf['treatment_specific_drug']))
    dic_primary_full, _ = return_metadata()

    dsiganalyzer = defaultdict(lambda: defaultdict(float))
    dsiganalyzer_simple = defaultdict(float)
    dsiganalyzer_simple_exposure = defaultdict(float)

    sigs_affected_signature_analyzer = signature_2_treatment_signature_analyzer[drug]
    for sample, d in df_sigan.iterrows():
        if sample in treated_samples[drug]['Pan']['YES']:
            dsiganalyzer[dic_primary_full[sample]][sample] = d.loc[sigs_affected_signature_analyzer].sum()
            dsiganalyzer_simple[sample] = d.loc[sigs_affected_signature_analyzer].sum()
            dsiganalyzer_simple_exposure[sample] = d.loc[sigs_affected_signature_analyzer].sum()/d.sum()

    dsigprof = defaultdict(lambda: defaultdict(float))
    dsigprof_simple = defaultdict(float)
    dsigprof_simple_exposure = defaultdict(float)
    sigs_affected_sigprofiler = signature_2_treatment_sigprofiler[drug]

    for sample, d in df_sigpro.iterrows():
        if sample in treated_samples[drug]['Pan']['YES']:
            all_treated_samples.add(sample)

            dsigprof[dic_primary_full[sample]][sample] = d.loc[sigs_affected_sigprofiler].sum()
            dsigprof_simple[sample] = d.loc[sigs_affected_sigprofiler].sum()
            dsigprof_simple_exposure[sample] = d.loc[sigs_affected_sigprofiler].sum()/d.sum()

    barplot_classes(dsigprof_simple_exposure, dsiganalyzer_simple_exposure, drug, outpath, 'count')
    plot_single_distribution(dsigprof_simple, dsiganalyzer_simple, drug, outpath, 'count')
    plot_single_correlation(dsigprof_simple, dsiganalyzer_simple, drug, outpath, 'count')
    plot_single_distribution(dsigprof_simple_exposure, dsiganalyzer_simple_exposure, drug, outpath, 'exposure')
    plot_single_correlation(dsigprof_simple_exposure, dsiganalyzer_simple_exposure, drug, outpath, 'exposure')

    bad_samples_drug = get_samples_similar_exposure(
        dsigprof_simple_exposure, dsiganalyzer_simple_exposure, drug, outpath, 'exposure'
    )


config_params(5)


def run():

    dic_primary_full, _ = return_metadata()
    signature_2_treatment_signature_analyzer = {
        'CISPLATIN': ['21_SBS31_0.953955_1', '14_1'],
        'CARBOPLATIN': ['21_SBS31_0.953955_1', '25_1'],
        '5-FU_CAPE': ['31_SBS17b_0.968799_1'],
        'OXALIPLATIN': ['14_1', '37_1']
    }
    signature_2_treatment_sigprofiler = {
        'CISPLATIN': ['1_SBS31_0.968153_0.98'],
        'CARBOPLATIN': ['1_SBS31_0.968153_0.98'],
        'OXALIPLATIN': ['20_0.92'],
        '5-FU_CAPE': ['19_SBS17b_0.961548_0.99'],
    }

    outpath = 'figures/'
    path_base = 'data/hartwig/signatures/extraction/results/{}/snvs/exposures/'
    path_siganalyzer = path_base.format("SignatureAnalyzer") + "Pan_full/Pan_full.exposures.tsv"
    path_sigprofiler = path_base.format("SigProfiler") + "PanNoSkinNoUnknown.snvs/PanNoSkinNoUnknown.snvs.exposures.tsv"

    for drug in signature_2_treatment_signature_analyzer:
        process_exposures_methods(
            path_siganalyzer, path_sigprofiler, drug,
            signature_2_treatment_signature_analyzer,
            signature_2_treatment_sigprofiler, outpath
        )


if __name__ == '__main__':
    run()
