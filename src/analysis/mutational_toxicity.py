# Import modules
import os
import gzip
import sys
from collections import defaultdict

import pandas as pd
import dill as pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import matplotlib as mpl

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from utils.colors_ttype import return_colors
from utils.metadata import return_metadata
from utils.plots_utils import config_params
from config import Conf

# Fix the random seed
np.random.seed(12345)

def return_sorted_specificmethod_dict(samples_mutations_treatment, df_muts, dic_primary_full):
    ttype_sorted = defaultdict(list)
    ttype_sorted_proportion = defaultdict(list)
    ttype_sorted_sample = defaultdict(list)
    ttype_sorted_proportion_sample = defaultdict(list)
    ttype_proportion_exposed = defaultdict(lambda: defaultdict(int))

    for sample, vals_mut in samples_mutations_treatment.items():
        # the total number of mutations is equivalent in both methods
        total_muts = df_muts.loc[sample].sum()
        ttype_sample = dic_primary_full[sample]

        if vals_mut > 0:
            ttype_sorted_proportion[ttype_sample].append(vals_mut / total_muts)
            ttype_sorted[ttype_sample].append(vals_mut)
            ttype_sorted_sample[sample] = vals_mut
            ttype_sorted_proportion_sample[sample] = vals_mut / total_muts
            ttype_proportion_exposed[ttype_sample]['exposed'] += 1

        else:
            ttype_proportion_exposed[ttype_sample]['not_exposed'] += 1

    sorted_ttypes = sorted(ttype_sorted, key=lambda k: np.median(ttype_sorted[k]), reverse=True)

    return sorted_ttypes, ttype_sorted, ttype_sorted_proportion, ttype_proportion_exposed


def return_sorted_coherent_dict(samples_to_analyze_coherent, df_sigprofiler, dic_primary_full, samples_mutations_treatment_sigan,
                                samples_mutations_treatment_sigpro):
    ttype_sorted = defaultdict(list)
    ttype_sorted_proportion = defaultdict(list)
    ttype_sorted_sample = defaultdict(list)
    ttype_sorted_proportion_sample = defaultdict(list)

    for sample in samples_to_analyze_coherent:

        # the total number of mutations is equivalent in both methods
        total_muts = df_sigprofiler.loc[sample].sum()
        ttype_sample = dic_primary_full[sample]
        mutation_toxicity_sigan = samples_mutations_treatment_sigan[sample]
        mutation_toxicity_sigpro = samples_mutations_treatment_sigpro[sample]

        mean_val = np.mean([mutation_toxicity_sigan, mutation_toxicity_sigpro])

        # plot only values higher than 0
        if mean_val > 0:
            ttype_sorted_proportion[ttype_sample].append(mean_val / total_muts)
            ttype_sorted[ttype_sample].append(mean_val)
            ttype_sorted_sample[sample] = mean_val
            ttype_sorted_proportion_sample[sample] = mean_val / total_muts

    sorted_ttypes = sorted(ttype_sorted, key=lambda k: np.median(ttype_sorted[k]), reverse=True)

    return sorted_ttypes, ttype_sorted, ttype_sorted_proportion, ttype_sorted_sample


def check_concordance_methods(proportion_treatment_mutations_sigprofiler, proportion_treatment_mutations_siganalyzer):
    if (proportion_treatment_mutations_sigprofiler == 0) & (proportion_treatment_mutations_siganalyzer < 0.0075):
        class_sample = 'not-exposed'

    elif (proportion_treatment_mutations_sigprofiler > 0) & (proportion_treatment_mutations_siganalyzer < 0.0075):
        class_sample = 'not-concordant'

    elif (proportion_treatment_mutations_sigprofiler == 0) & (proportion_treatment_mutations_siganalyzer >= 0.0075):
        class_sample = 'not-concordant'

    elif (proportion_treatment_mutations_sigprofiler > 0) & (proportion_treatment_mutations_siganalyzer >= 0.0075):

        if abs(proportion_treatment_mutations_sigprofiler - proportion_treatment_mutations_siganalyzer) <= 0.15:
            class_sample = 'concordant'
        else:
            class_sample = 'not-concordant'

    return class_sample

def plot_sigmoid_tumortype(wanted_ttypes, sorted_ttypes, ttype_sorted, outname='full'):

    color_ttype = return_colors()
    # now plot the distribution of mutations
    config_params(5)

    fig, axs = plt.subplots(1, len(wanted_ttypes), figsize=(1.75, 1.5), sharey=True)  # ,  sharey='row')
    plt.yscale('log')

    toplot = 0

    for ttype in sorted_ttypes:

        if ttype in wanted_ttypes:
            print('{} ({})'.format(ttype, len(ttype_sorted[ttype])))
            axs[toplot].set_xlabel('{} ({})'.format(ttype, len(ttype_sorted[ttype])), rotation=90)
            axs[toplot].set_ylim(10, 200000)
            range_nums = [num / len(ttype_sorted[ttype]) for num in range(len(ttype_sorted[ttype]))]
            axs[toplot].hlines(np.median(ttype_sorted[ttype]),
                               np.median(range_nums) - 0.25, np.median(range_nums) + 0.45,
                               color='grey', alpha=0.4)
            axs[toplot].scatter(range_nums,
                                sorted(ttype_sorted[ttype]), s=2,
                                linewidth=0.2, color=color_ttype[ttype],
                                )

            axs[toplot].set_xlim(-0.2, 1.1)

            axs[toplot].spines['bottom'].set_visible(False)
            axs[toplot].spines['right'].set_visible(False)
            axs[toplot].set_xticklabels([])

            if toplot > 0:
                axs[toplot].spines['left'].set_visible(False)
            else:
                axs[toplot].set_ylabel('Number of SBS\nassociated to treatments')

            axs[toplot].xaxis.set_label_position('top')
            axs[toplot].xaxis.set_ticks_position('none')
            axs[toplot].yaxis.set_ticks_position('none')

            toplot += 1
    plt.savefig('figures/fig5_sigmoid_{}.svg'.format(outname))
    plt.close()


def plot_sigmoid_proportion(wanted_ttypes, sorted_ttypes, ttype_sorted_proportion, outname='full'):

    config_params(5)
    color_ttype = return_colors()

    fig, axs = plt.subplots(1, len(wanted_ttypes), figsize=(1.75, 1.5), sharey=True)

    toplot = 0
    for ttype in sorted_ttypes:

        if ttype in wanted_ttypes:

            range_nums = [num / len(ttype_sorted_proportion[ttype]) for num in range(len(ttype_sorted_proportion[ttype]))]

            axs[toplot].scatter(range_nums,
                                sorted(ttype_sorted_proportion[ttype]), s=3,
                                linewidth=0.2, color=color_ttype[ttype],
                                )

            axs[toplot].hlines(np.median(ttype_sorted_proportion[ttype]), np.median(range_nums) - 0.25, np.median(range_nums) + 0.45,
                               color='grey', alpha=0.4)  # , s = 10)

            axs[toplot].spines['top'].set_visible(False)
            axs[toplot].spines['right'].set_visible(False)
            axs[toplot].set_xlim(-0.2, 1.5)

            axs[toplot].set_xticklabels([])

            if toplot > 0:

                axs[toplot].spines['left'].set_visible(False)
            else:
                axs[toplot].set_ylabel('Proportion of SBS\nassociated to treatments')
                axs[toplot].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
                axs[toplot].set_yticklabels(["0", "20", "40", "60", "80", "100"])

            axs[toplot].xaxis.set_ticks_position('none')
            axs[toplot].yaxis.set_ticks_position('none')

            toplot += 1

    plt.ylim(0, 0.85)
    plt.savefig('figures/fig5_proportion_{}.svg'.format(outname))
    plt.close()


def mutational_toxicity_per_sample(df_siganalyzer, df_sigprofiler, siganalyzer_sigs, sigprofiler_sigs, total_samples, list_samples_sigprofiler,
                                   treated_samples_dictionary):
    # this keeps track on whether the signature has already been assessed for the sample
    sample_sigs_examined_sigpro = defaultdict(set)
    sample_sigs_examined_sigan = defaultdict(set)

    # mutational burden or toxicity per method
    samples_mutations_treatment_sigpro = defaultdict(float)
    samples_mutations_treatment_sigan = defaultdict(float)

    # now we should assess the overall mutation burden, removing the forbidden and taking
    # into account only once the sigatures related to the treatments

    for drug in sigprofiler_sigs:

        if drug != 'AGING':

            samples_to_analyze_sigprofiler = [s for s in treated_samples_dictionary[drug]['Pan']['YES'] if s in list_samples_sigprofiler]
            samples_to_analyze_siganalyzer = [s for s in treated_samples_dictionary[drug]['Pan']['YES'] if s in total_samples]

            for sample in samples_to_analyze_siganalyzer:

                sample_exposure = df_siganalyzer.loc[sample]
                all_signatures = siganalyzer_sigs[drug]

                for sig in all_signatures:
                    # this checks that the sample does not have a repeated count, because it might have been under two treatments with same sig
                    if sig not in sample_sigs_examined_sigan[sample]:
                        samples_mutations_treatment_sigan[sample] += round(sample_exposure[sig].sum())
                        sample_sigs_examined_sigan[sample].add(sig)

            for sample in samples_to_analyze_sigprofiler:

                sample_exposure = df_sigprofiler.loc[sample]
                all_signatures = sigprofiler_sigs[drug]

                for sig in all_signatures:

                    # this checks that the sample does not have a repeated count
                    if sig not in sample_sigs_examined_sigpro[sample]:
                        samples_mutations_treatment_sigpro[sample] += round(sample_exposure[sig].sum())
                        sample_sigs_examined_sigpro[sample].add(sig)

    return samples_mutations_treatment_sigpro, samples_mutations_treatment_sigan


# Fig5b and EDF9b, 10b
def sigmoid_plots(wanted_ttypes, allowed_samples, forbidden_samples, dic_primary_full, df_siganalyzer, df_sigprofiler, siganalyzer_sigs,
                  sigprofiler_sigs, total_samples, list_samples_sigprofiler, treated_samples_dictionary):

    # -----
    # Figure 5b and Extended8d and 9d
    # -----
    samples_mutations_treatment_sigpro, samples_mutations_treatment_sigan = mutational_toxicity_per_sample(df_siganalyzer, df_sigprofiler,
                                                                                                           siganalyzer_sigs, sigprofiler_sigs,
                                                                                                           total_samples, list_samples_sigprofiler,
                                                                                                           treated_samples_dictionary)

    # the samples to analyze in thecCoherent scenario
    samples_to_analyze_coherent = [s for s in allowed_samples if s not in forbidden_samples]

    # Coherent plot plot
    sorted_ttypes, ttype_sorted, ttype_sorted_proportion, samples_mutations_treatment_coherent = return_sorted_coherent_dict(
        samples_to_analyze_coherent, df_sigprofiler, dic_primary_full, samples_mutations_treatment_sigan, samples_mutations_treatment_sigpro)

    plot_sigmoid_tumortype(wanted_ttypes, sorted_ttypes, ttype_sorted, 'coherent')
    plot_sigmoid_proportion(wanted_ttypes, sorted_ttypes, ttype_sorted_proportion, 'coherent')

    # SignatureAnalyzer plot
    sorted_ttypes, ttype_sorted, ttype_sorted_proportion, ttype_proportion_exposed = return_sorted_specificmethod_dict(
        samples_mutations_treatment_sigan,
        df_siganalyzer, dic_primary_full)

    sorted_wanted_ttypes = [s for s in sorted_ttypes if s in wanted_ttypes]
    plot_sigmoid_tumortype(wanted_ttypes, sorted_ttypes, ttype_sorted, 'SignatureAnalyzer')
    plot_sigmoid_proportion(wanted_ttypes, sorted_ttypes, ttype_sorted_proportion, 'SignatureAnalyzer')
    plot_barplot_proportion_method(ttype_proportion_exposed, sorted_wanted_ttypes, 'SignatureAnalyzer')
    # SigProfiler plot
    sorted_ttypes, ttype_sorted, ttype_sorted_proportion, ttype_proportion_exposed = return_sorted_specificmethod_dict(
        samples_mutations_treatment_sigpro,
        df_sigprofiler, dic_primary_full)

    sorted_wanted_ttypes = [s for s in sorted_ttypes if s in wanted_ttypes]
    plot_sigmoid_tumortype(wanted_ttypes, sorted_ttypes, ttype_sorted, 'SigProfiler')
    plot_sigmoid_proportion(wanted_ttypes, sorted_ttypes, ttype_sorted_proportion, 'SigProfiler')
    plot_barplot_proportion_method(ttype_proportion_exposed, sorted_wanted_ttypes, 'SigProfiler')

    return samples_mutations_treatment_coherent, samples_mutations_treatment_sigpro, samples_mutations_treatment_sigan


def return_IQR(x):
    q75, q25 = np.percentile(x, [75, 25])
    iqr = q75 - q25

    return iqr


def color_treatments_f():
    color_treatments = {
        '5-FU_CAPE': '#2c89a0ff',
        'CAPECITABINE': '#2c89a0ff',
        'OXALIPLATIN': '#ffb380ff',
        'CISPLATIN': '#ffe6d5ff',
        'CARBOPLATIN': '#ff6600ff',
        'AGING': 'green',
        'SMOKING': 'black'
    }
    return color_treatments


def plot_wiskeys(wanted_ttypes, coherent_samples_ttype, type_extraction):

    config_params(5.5)

    color_treatments = color_treatments_f()

    fig, axs = plt.subplots(1, len(wanted_ttypes), figsize=(6, 0.8))
    toplot = 0

    for ttype in wanted_ttypes:
        labels = []
        maxval = []
        count_pos = 0
        for ix, drug in enumerate(['AGING', 'CAPECITABINE', '5-FU_CAPE', 'CISPLATIN', 'CARBOPLATIN', 'OXALIPLATIN']):
            ttype_dict = coherent_samples_ttype[drug][ttype]
            if len(ttype_dict['count']) > 10:
                median_val = np.median(ttype_dict['count'])
                if median_val > 0:
                    iqr_range = return_IQR(ttype_dict['count'])
                    maxval.append(iqr_range + median_val)
                    axs[toplot].vlines(count_pos, median_val - iqr_range, median_val + iqr_range, lw=1,
                                       color=color_treatments[drug])

                axs[toplot].scatter(count_pos, median_val, color=color_treatments[drug], label=drug,
                                    linewidths=1)
                labels.append(len(ttype_dict['count']))
                count_pos += 1

        axs[toplot].set_title(ttype)
        axs[toplot].set_xticks([i for i in range(len(labels))])
        axs[toplot].set_xticklabels(labels, ha='center')
        axs[toplot].set_xlim(-1, len(labels))
        axs[toplot].set_ylim(0, np.max(maxval))
        axs[toplot].spines['top'].set_visible(False)
        axs[toplot].spines['right'].set_visible(False)
        toplot += 1

    plt.legend()
    os.makedirs('figures/{}'.format(type_extraction), exist_ok=True)
    plt.savefig('figures/{}/wiskey.svg'.format(type_extraction))
    plt.close()


def plot_barplots(wanted_ttypes, dexposures_class):
    color_treatments = color_treatments_f()
    config_params(4.5)
    fig, axs = plt.subplots(1, len(wanted_ttypes), figsize=(4.9, 0.35))

    toplot = 0
    for ttype in wanted_ttypes:

        labels = []
        count_pos = 0

        for ix, treat in enumerate(['AGING', 'CAPECITABINE', '5-FU_CAPE', 'CISPLATIN', 'CARBOPLATIN', 'OXALIPLATIN']):
            d_class = dexposures_class[treat][ttype]
            if sum(list(d_class.values())) > 10:
                concordant = d_class.get('concordant', 0)
                not_concordant = d_class.get('not-concordant', 0)
                not_exposed = d_class.get('not-exposed', 0)
                total = concordant + not_concordant + not_exposed

                if concordant > 0:
                    Nconcordant = concordant / total
                    Nnot_concordant = not_concordant / total
                    Nnot_exposed = not_exposed / total
                    axs[toplot].bar(count_pos, Nconcordant, color=color_treatments[treat],
                                    linewidth=0.5, edgecolor='black', )
                    axs[toplot].bar(count_pos, Nnot_concordant, bottom=Nconcordant, color=color_treatments[treat],
                                    linewidth=0.5, edgecolor='black', alpha=0.3,
                                    hatch="//////")
                    axs[toplot].bar(count_pos, Nnot_exposed, bottom=(Nconcordant + Nnot_concordant), color='#f0f0f0',
                                    linewidth=0.5, edgecolor='black', alpha=0.3, )
                    count_pos += 1

                    labels.append(total)

        axs[toplot].set_xticks([i for i in range(len(labels))])
        axs[toplot].set_xticklabels(labels, ha='center')
        axs[toplot].set_xlim(-1, len(labels))
        axs[toplot].spines['top'].set_visible(False)
        axs[toplot].set_ylim(0, 1)
        axs[toplot].spines['right'].set_visible(False)
        toplot += 1

    plt.savefig('figures/wikeys_coherence_box.svg')

    plt.close()


def plot_barplot_proportion_method(ttype_proportion_exposed, sorted_wanted_ttypes, type_extraction):

    config_params(5)

    fig, axs = plt.subplots(1, len(sorted_wanted_ttypes), figsize=(2.5, 0.5), sharey=True)

    col_ttype = return_colors()
    toplot = 0

    for ttype in sorted_wanted_ttypes:

        exp_ttype = ttype_proportion_exposed[ttype]['exposed']
        notexp_ttype = ttype_proportion_exposed[ttype]['not_exposed']
        axs[toplot].set_xlabel('{}'.format(exp_ttype + notexp_ttype))

        perc = 100 * exp_ttype / (exp_ttype + notexp_ttype)
        axs[toplot].bar(toplot, perc, color=col_ttype[ttype],
                        )

        axs[toplot].bar(toplot, 100 - perc, bottom=perc, color='#f0f0f0',
                        )
        axs[toplot].xaxis.set_ticks_position('none')
        axs[toplot].yaxis.set_ticks_position('none')
        axs[toplot].text(toplot - 0.20, 110, int(100 * round(exp_ttype / (exp_ttype + notexp_ttype), 2)))

        axs[toplot].spines['top'].set_visible(False)
        # ax.spines['left'].set_visible(False)
        if toplot > 0:
            axs[toplot].spines['left'].set_visible(False)

        axs[toplot].spines['right'].set_visible(False)
        axs[toplot].set_xticklabels([])

        toplot += 1
    os.makedirs('figures/{}'.format(type_extraction), exist_ok=True)
    plt.savefig('figures/{}/supp_barplot.svg'.format(type_extraction))
    plt.close()


def return_aging():
    old_metadata = Conf['age_metadata']

    meta = pd.read_csv(old_metadata, sep='\t')

    meta['yearbiopsy'] = meta['biopsyDate'].apply(lambda x: float(str(x).split('-')[0]))

    meta['AGE'] = meta['yearbiopsy'] - meta['birthYear']

    dic_age = dict(zip(meta.sampleId, meta.AGE))

    return dic_age


def proportion_aging_signature(dexposures_specific_sample, counter_length, dic_age, dic_primary_full):
    dic_results = defaultdict(lambda: defaultdict(dict))
    dic_results_ttype = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    color_ttype = return_colors()

    for drug in dexposures_specific_sample.keys():

        sig1 = []
        treat = []
        total_1 = []
        total_treat = []
        quocient = []
        colors_1 = []

        if drug != 'AGING':

            for sample, nsigdrug in dexposures_specific_sample[drug].items():

                # check if we also have data on the drug we are looking at
                if drug in counter_length[sample].keys():

                    # check the total days of treatment
                    total_days_treatment = [i for i in counter_length[sample][drug] if str(i) != 'nan']

                    # sum all days of exposure
                    sum_days = np.sum(total_days_treatment)

                    if sum_days > 0:

                        # check if the sample has age info
                        if sample in dic_age.keys():
                            years_life = dic_age[sample]

                            # sanity check
                            if years_life > 0:

                                nsig1 = dexposures_specific_sample['AGING'].get(sample, 0)
                                if nsig1 > 0:

                                    months_life = years_life * 12
                                    norm_sig1 = nsig1 / months_life

                                    if (nsig1 > 1) & (nsigdrug > 1):

                                        months_treat = sum_days / 31
                                        norm_sigdrug = nsigdrug / months_treat

                                        sig1.append(norm_sig1)

                                        treat.append(norm_sigdrug)
                                        total_1.append(nsig1)
                                        total_treat.append(nsigdrug)
                                        ttype = dic_primary_full[sample]

                                        dic_results_ttype['AGING'][ttype][sample] = norm_sig1
                                        dic_results_ttype[drug][ttype]['normtreat'].append(norm_sigdrug)
                                        dic_results_ttype[drug][ttype]['totalsig1'].append(nsig1)
                                        dic_results_ttype[drug][ttype]['totaltreat'].append(nsigdrug)
                                        dic_results_ttype[drug][ttype]['treated_days'].append(sum_days)

                                        colors_1.append(color_ttype[ttype])

                                        if norm_sig1 > 0:
                                            quocient.append(norm_sigdrug / norm_sig1)

            dic_results[drug]['normsig1'] = sig1
            dic_results[drug]['normtreat'] = treat
            dic_results[drug]['totalsig1'] = total_1
            dic_results[drug]['totaltreat'] = total_treat
            dic_results[drug]['colors'] = colors_1

    return dic_results, dic_results_ttype


def plot_aging_twosigs(dic_results, method_extraction):
    for drug in dic_results:
        sig1 = dic_results[drug]['normsig1']
        treat = dic_results[drug]['normtreat']
        total_1 = dic_results[drug]['totalsig1']
        total_treat = dic_results[drug]['totaltreat']
        colors = dic_results[drug]['colors']

        config_params(5.5)

        fig, axs = plt.subplots(1, 2, figsize=(1.7, 1))
        g = sns.boxplot(
            data=[sig1, treat], palette=['darkred', '#cbcacb'],
            ax=axs[1], showfliers=False, linewidth=0.5
        )
        g.set_yscale('log')

        axs[1].scatter(
            [0 + np.random.uniform(-0.2, 0.2, 1)[0] for i in sig1],
            sig1, color=colors, s=0.5, alpha=0.4
        )
        axs[1].scatter(
            [1 + np.random.uniform(-0.2, 0.2, 1)[0] for i in treat],
            treat, color=colors, s=0.5, alpha=0.4
        )

        axs[1].set_ylabel('Mutations per\nmonth of exposure', fontsize=6)
        axs[1].set_xticklabels(
            ['{}'.format(len(total_1)),
             '{}'.format(len(total_treat))], rotation=90, fontsize=6
        )
        axs[1].spines['right'].set_visible(False)
        axs[1].spines['top'].set_visible(False)

        g = sns.boxplot(data=[total_1, total_treat], palette=['darkred', '#cbcacb'], ax=axs[0],
                        showfliers=False, linewidth=0.5)
        g.set_yscale('log')
        axs[0].scatter(
            [0 + np.random.uniform(-0.2, 0.2, 1)[0] for i in total_1],
            total_1, color=colors, s=0.5, alpha=0.4
        )
        axs[0].scatter(
            [1 + np.random.uniform(-0.2, 0.2, 1)[0] for i in total_treat],
            total_treat, color=colors, s=0.5, alpha=0.4
        )
        axs[0].set_xticklabels([])
        axs[0].spines['right'].set_visible(False)
        axs[0].spines['top'].set_visible(False)
        axs[0].set_ylabel('Mutations', fontsize=6)
        plt.suptitle(drug)
        plt.savefig('figures/boxplot_months_of_exposure_{}.svg'.format(method_extraction))
        plt.close()


def plot_lines_probabilities(prob_ttypes, method_extraction):
    color_treatments = color_treatments_f()

    wanted_ttypes_2 = ['Colon-Rectum', 'Breast', 'Lung', 'Ovary']

    for ttype in wanted_ttypes_2:

        fig, ax = plt.subplots(1, 1, figsize=(1.1, 0.8))
        config_params(5.5)

        plt.title(ttype)
        for treatment in prob_ttypes[ttype]:
            plt.plot(
                [0, 2500], [0, 2500 * prob_ttypes[ttype][treatment]['gene']],
                lw=2, label=treatment, color=color_treatments[treatment]
            )
            if treatment == 'AGING':
                plt.hlines(
                    1000 * prob_ttypes[ttype][treatment]['gene'],
                    0, 1000, lw=1, linestyle='--', color='grey'
                )
                plt.vlines(
                    1000, 0, 1000 * prob_ttypes[ttype][treatment]['gene'],
                    lw=1, linestyle='--', color='grey'
                )

            print(ttype, treatment, 1000 * prob_ttypes[ttype][treatment]['gene'],
                  1000 * prob_ttypes[ttype][treatment]['cgc'], )

        # plt.legend(bbox_to_anchor = [1.1, 1.1])
        plt.ylabel('Risk of coding affecting')
        plt.xlabel('Mutations')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.savefig('figures/lines_probabilities_{}.svg'.format(method_extraction))
        plt.close()


def plot_lines_time(prob_ttypes, dic_results_ttype, coding_type, method_extraction):
    color_treatments = color_treatments_f()

    ttypes_selected = ['Colon-Rectum', 'Breast', 'Ovary']
    total = []

    for ttype in prob_ttypes:

        for treatment in prob_ttypes[ttype].keys():
            treatment_fixed = treatment

            if treatment == 'CAPECITABINE':
                treatment_fixed = '5-FU_CAPE'

            if treatment_fixed not in ['AGING', 'SMOKING']:

                for sample in dic_results_ttype[treatment_fixed][ttype]['treated_days']:
                    total.append(sample / 7)

    ttypes_selected.append('Lung')

    for ttype in ttypes_selected:

        fig, axs = plt.subplots(
            nrows=1, ncols=2, figsize=(2.75, 1.5), sharey=True,
            gridspec_kw={'width_ratios': [3, 2]}
        )
        config_params(5.5)
        axs[0].set_xticks(np.arange(0, 52, 8))

        plt.yscale('log')
        toplot_aging = []

        for treatment in prob_ttypes[ttype].keys():

            treatment_fixed = treatment

            if treatment == 'CAPECITABINE':
                treatment_fixed = '5-FU_CAPE'

            if treatment not in ['AGING', 'SMOKING']:
                prob_treatment = prob_ttypes[ttype][treatment][coding_type]

                if len(dic_results_ttype[treatment_fixed][ttype]['normtreat']) > 1:
                    median_val = np.median(dic_results_ttype[treatment_fixed][ttype]['treated_days']) / 7

                    if median_val > 0:
                        toplot_aging.append(median_val)

                        print(treatment, ttype, len(dic_results_ttype[treatment_fixed][ttype]['normtreat']))
                        axs[0].vlines(median_val, 0, 80, lw=0.6, color='grey', linestyles='dashed')
                        sig_treat_high = 7 * (np.percentile(
                            dic_results_ttype[treatment_fixed][ttype]['normtreat'], 75)
                        ) / 31
                        sig_treat_low = 7 * (np.percentile(
                            dic_results_ttype[treatment_fixed][ttype]['normtreat'], 25)
                        ) / 31
                        sig_treat_median = 7 * (np.median(
                            dic_results_ttype[treatment_fixed][ttype]['normtreat'])
                        ) / 31

                        weeks_to_plot = 40

                        axs[0].fill_between(
                            x=np.arange(weeks_to_plot),
                            y1=[prob_treatment * i * sig_treat_high for i in np.arange(weeks_to_plot)],
                            y2=[prob_treatment * i * sig_treat_low for i in np.arange(weeks_to_plot)],
                            color=color_treatments[treatment_fixed],
                            alpha=0.3,
                            label=treatment_fixed
                        )

                        axs[0].scatter(
                            x=median_val,
                            y=prob_treatment * median_val * sig_treat_median,
                            s=20, c=color_treatments[treatment_fixed],
                            alpha=1,
                            linewidths=0.5,
                            edgecolors='black'
                        )

        treatment = 'AGING'

        prob_aging = prob_ttypes[ttype]['AGING'][coding_type]
        sig1_per_mont = np.array(list(dic_results_ttype['AGING'][ttype].values()))
        print(ttype, 'SIG1', len(sig1_per_mont))
        sig1_high = 7 * (np.percentile(sig1_per_mont, 75)) / 30
        sig1_low = 7 * (np.percentile(sig1_per_mont, 25)) / 30
        sig_treat_median = 7 * (np.median(sig1_per_mont)) / 31

        weeks_to_plot = 40

        axs[0].fill_between(
            x=np.arange(weeks_to_plot),
            y1=[prob_aging * i * sig1_high for i in np.arange(weeks_to_plot)],
            y2=[prob_aging * i * sig1_low for i in np.arange(weeks_to_plot)],
            color=color_treatments[treatment],
            alpha=0.3,
            label=treatment
        )

        for med in toplot_aging:
            axs[0].scatter(
                x=med,
                y=prob_aging * med * sig_treat_median,
                s=20, c=color_treatments[treatment],
                alpha=1,
                linewidths=0.5,
                edgecolors='black'
            )
        years_to_plot = 11

        sig1_high = 12 * (np.percentile(sig1_per_mont, 75))
        sig1_low = 12 * (np.percentile(sig1_per_mont, 25))
        sig_treat_median = 12 * (np.median(sig1_per_mont))

        axs[1].fill_between(
            x=np.arange(1, years_to_plot),
            y1=[prob_aging * i * sig1_high for i in np.arange(1, years_to_plot)],
            y2=[prob_aging * i * sig1_low for i in np.arange(1, years_to_plot)],
            color=color_treatments[treatment],
            alpha=0.3,
            label=treatment
        )

        axs[1].scatter(
            x=np.arange(1, years_to_plot, 4),
            y=[prob_aging * i * sig_treat_median for i in np.arange(1, years_to_plot, 4)],
            s=10,
            c=color_treatments[treatment],
            alpha=1
        )

        axs[0].set_ylabel('Risk of {} affecting'.format(coding_type))
        axs[0].set_xlabel('Weeks of exposure')
        axs[1].set_xticks(np.arange(1, years_to_plot, 1))
        axs[1].set_xlim(0.5, years_to_plot)
        axs[1].set_xlabel('Years of exposure')
        axs[1].spines['top'].set_visible(False)
        axs[1].spines['left'].set_visible(False)
        axs[0].spines['top'].set_visible(False)
        axs[0].spines['right'].set_visible(False)
        axs[1].yaxis.set_ticks_position('none')
        axs[1].spines['right'].set_visible(False)
        axs[0].set_yticks([0.01, 0.1, 1, 10, 100])
        axs[0].set_yticklabels(["0.01", "0.1", "1", "10", "100"])

        plt.title(ttype)
        plt.tight_layout()
        axs[0].legend()

        plt.savefig('figures/probability_landscape_{}_{}.svg'.format(method_extraction, coding_type))
        plt.close()


def treatments_literature():
    dic_treatments_norm_time = {
        'Colon-Rectum': {
            'Capecitabine/Oxaliplatin': 3,
            'Bevacizumab/Capecitabine/Oxaliplatin': 3,
            'Bevacizumab/Capecitabine/Capecitabine/Oxaliplatin': 3,
            'Bevacizumab/Capecitabine/Irinotecan/Oxaliplatin': 3,
            'Capecitabine/Irinotecan/Oxaliplatin': 3,
            'Bevacizumab/Capecitabine': 3,
            'Capecitabine/Mitomycin C': 3,
            'Bevacizumab/Capecitabine/Capecitabine/Irinotecan/Oxaliplatin': 3
        },
        'Ovary': {
            'Carboplatin/Paclitaxel': 3,
            'Carboplatin/Carboplatin/Paclitaxel/Paclitaxel': 3,
            'Bevacizumab/Carboplatin/Carboplatin/Gemcitabine/Paclitaxel': 3,

        },
        'Lung': {
            'Cisplatin/Pemetrexed': 3,
            'Carboplatin/Pemetrexed': 3,
            'Cisplatin': 3,
            'Cisplatin/Etoposide': 3,
            'Cisplatin/Nivolumab/Pemetrexed': 4,
            'Carboplatin/Gemcitabine': 3,
        },
        'Breast': {
            'Cyclophosphamide/Docetaxel/Epirubicin/Fluorouracil/Tamoxifen': 3,
            'Anastrozole/Cyclophosphamide/Epirubicin/Fluorouracil/Tamoxifen': 3,
            'Capecitabine': 6,
            'Carboplatin': 3
        }
    }

    return dic_treatments_norm_time


def select_treated_samples_regime(df, dic_treatments_norm_time):

    if df['primaryTumorLocation_fixed'] == 'Breast':
        treat = str(df['preTreatments'])
        if 'Capecitabine' in treat:
            time_to_return = 6
        elif ('Carboplatin' in treat) or ('Fluorouracil' in treat):
            time_to_return = 3
        else:
            time_to_return = np.nan
    else:
        time_to_return = dic_treatments_norm_time[df['primaryTumorLocation_fixed']].get(
            df['preTreatments_fixed'], np.nan
        )

    return time_to_return


def get_complete_regime():
    dic_treatments_norm_time = treatments_literature()
    metadata = pd.read_csv('data/clinical_data/preTreatments_fixed.tsv', sep='\t')

    # we will get only the samples with the specific combination of treatments, except for Breast
    metadata_small = metadata[metadata['primaryTumorLocation_fixed'].isin(dic_treatments_norm_time)]
    metadata_small['time_treat'] = metadata_small.apply(
        select_treated_samples_regime, args=(dic_treatments_norm_time,), axis=1
    )
    metadata_small = metadata_small[~metadata_small['time_treat'].isnull()]

    dict_treatments_literature = dict(zip(metadata_small['sampleId'], metadata_small['preTreatments_fixed']))
    dict_ttype = dict(zip(metadata_small['sampleId'], metadata_small['primaryTumorLocation_fixed']))
    dict_treatments_literature_time = dict(zip(metadata_small['sampleId'], metadata_small['time_treat']))

    dic_lit = defaultdict(lambda: defaultdict(int))

    dic_equi = {'Fluorouracil': '5-FU_CAPE', 'Capecitabine': '5-FU_CAPE'}

    for sample, treatment in dict_treatments_literature.items():
        for treat in treatment.split('/'):
            treat = dic_equi.get(treat, treat)
            if (dict_ttype[sample] == 'Breast') & (treat in ['Capecitabine', 'Fluorouracil']):
                dic_lit[sample]['5-FU_CAPE'] = [dict_treatments_literature_time[sample] * 31]
            elif (dict_ttype[sample] == 'Breast') & (treat == 'Carboplatin'):
                dic_lit[sample]['CARBOPLATIN'] = [dict_treatments_literature_time[sample] * 31]
            else:
                dic_lit[sample][treat.upper()] = [dict_treatments_literature_time[sample] * 31]

    return dic_lit


def write_supptable(df_siganalyzer, wanted_ttypes, forbidden_samples,
                    coherent_samples_samplelevel, exposures_drug_sigprofiler,
                    exposures_drug_siganalyzer, samples_mutations_treatment_coherent,
                    samples_mutations_treatment_sigpro, samples_mutations_treatment_sigan, dic_primary_full):
    # get index
    equivalent_ids = {v: ix for ix, v in enumerate(df_siganalyzer.index)}
    muts_per_sample = df_siganalyzer.sum(axis=1).to_dict()

    order_cols = [
        'Sample_ID', 'Tumor type', 'Proportion treatment related mutations', 'Total treatment related mutations',
        'AGING', '5-FU_CAPE', 'CARBOPLATIN', 'CISPLATIN', 'OXALIPLATIN',
        'Gene-affecting AGING', 'Gene-affecting 5-FU_CAPE', 'Gene-affecting CARBOPLATIN', 'Gene-affecting CISPLATIN',
        'Gene-affecting OXALIPLATIN',
        'CGC-affecting AGING', 'CGC-affecting 5-FU_CAPE', 'CGC-affecting CARBOPLATIN', 'CGC-affecting CISPLATIN',
        'CGC-affecting OXALIPLATIN',
    ]

    prob_ttypes = pickle.load(gzip.open('data/probability_MB/results_probability_sigpro_20190502.pckl.gz'))
    all_dics = [coherent_samples_samplelevel, exposures_drug_sigprofiler, exposures_drug_siganalyzer]
    all_treatment_dics = [
        samples_mutations_treatment_coherent, samples_mutations_treatment_sigpro, samples_mutations_treatment_sigan
    ]
    labels = ['Coherent', 'SigProfiler', 'SignatureAnalyzer']

    for ix, dic_expo in enumerate(all_dics):

        treatment_muts = all_treatment_dics[ix]
        # create table
        results_to_table = defaultdict(lambda: defaultdict(float))

        for treatment, d in dic_expo.items():

            treatm = treatment
            if treatment == '5-FU_CAPE':
                treatm = 'CAPECITABINE'
            for sample, values in d.items():
                if (ix != 0) or (sample not in forbidden_samples):
                    if dic_primary_full[sample] in wanted_ttypes:
                        ttype = dic_primary_full[sample]
                        results_to_table[sample][treatment] = values
                        results_to_table[sample]['Tumor type'] = ttype
                        prob_genes = prob_ttypes[ttype][treatm]['gene'] * values
                        prob_cgc = prob_ttypes[ttype][treatm]['cgc'] * values
                        results_to_table[sample]['Gene-affecting {}'.format(treatment)] = prob_genes
                        results_to_table[sample]['CGC-affecting {}'.format(treatment)] = prob_cgc
                        results_to_table[sample]['Total treatment related mutations'] = treatment_muts.get(
                            sample, np.nan
                        )
                        results_to_table[sample]['Proportion treatment related mutations'] = treatment_muts.get(
                            sample, np.nan
                        ) / muts_per_sample.get(sample, np.nan)
                        results_to_table[sample]['Sample_ID'] = equivalent_ids[sample]

        suppA = pd.DataFrame(results_to_table).T[order_cols]
        suppA.to_csv('figures/table_sheet_{}.tsv'.format(labels[ix]), sep='\t',
                     index=False, header=True, na_rep='NA')


def read_exposures(exposures_path):
    df = pd.read_csv(exposures_path, sep='\t', index_col=0)
    df = df.T
    return df


def analysis_coherence(
        df_siganalyzer, df_sigprofiler, siganalyzer_sigs,
        sigprofiler_sigs, treated_samples_dictionary, dic_primary_full):

    # we will annotate the samples that have non coherent exposures in the treatment signatures
    forbidden_samples = set()
    allowed_samples = set()

    class_samples_according_coherence = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    coherent_samples_ttype = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    coherent_samples_samplelevel = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    exposures_drug_sigprofiler = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    exposures_drug_siganalyzer = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    exposures_drug_sigprofiler_ttype = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    exposures_drug_siganalyzer_ttype = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    total_samples = df_siganalyzer.index.tolist()
    list_samples_sigprofiler = set(df_sigprofiler.index.tolist())

    # for each of the drug related to a signature
    for drug in tqdm(sigprofiler_sigs):

        sigs_sigan = siganalyzer_sigs[drug]
        sigs_sigpro = sigprofiler_sigs[drug]

        # SigProfiler analyzes less samples than SignatureAnalyzer
        for sample in total_samples:

            # if sample is treated, then get the value
            if (drug == 'AGING') or (sample in treated_samples_dictionary[drug]['Pan']['YES']):

                treatment_mutation_count_siganalyzer = df_siganalyzer.loc[sample][sigs_sigan].sum()

                proportion_treatment_mutations_siganalyzer = treatment_mutation_count_siganalyzer / df_siganalyzer.loc[sample].sum()
                exposures_drug_siganalyzer[drug][sample] = treatment_mutation_count_siganalyzer

                if proportion_treatment_mutations_siganalyzer >= 0.0075:
                    exposures_drug_siganalyzer_ttype[drug][dic_primary_full[sample]]['count'].append(
                        treatment_mutation_count_siganalyzer
                    )
                    exposures_drug_siganalyzer_ttype[drug][dic_primary_full[sample]]['exposure'].append(
                        proportion_treatment_mutations_siganalyzer
                    )

                treatment_mutation_count_sigprofiler = None
                proportion_treatment_mutations_sigprofiler = None

                if sample in list_samples_sigprofiler:
                    treatment_mutation_count_sigprofiler = df_sigprofiler.loc[sample][sigs_sigpro].sum()
                    proportion_treatment_mutations_sigprofiler = treatment_mutation_count_sigprofiler / df_sigprofiler.loc[sample].sum()
                    exposures_drug_sigprofiler[drug][sample] = treatment_mutation_count_sigprofiler

                    if treatment_mutation_count_sigprofiler > 0:
                        exposures_drug_sigprofiler_ttype[drug][dic_primary_full[sample]]['count'].append(
                            treatment_mutation_count_sigprofiler
                        )
                        exposures_drug_sigprofiler_ttype[drug][dic_primary_full[sample]]['exposure'].append(
                            proportion_treatment_mutations_sigprofiler
                        )

                # if the sample is also available in SigProfiler
                if proportion_treatment_mutations_sigprofiler is not None:

                    class_sample = check_concordance_methods(
                        proportion_treatment_mutations_sigprofiler,
                        proportion_treatment_mutations_siganalyzer
                    )
                    # annotate the class of coherence in this particular drug
                    class_samples_according_coherence[drug][dic_primary_full[sample]][class_sample] += 1

                    if class_sample == 'concordant':
                        if drug != 'AGING':
                            # this means we have data on exposures for both methods and
                            # in a specific case it was concordant
                            allowed_samples.add(sample)

                        mean_count_methods = np.mean([
                            treatment_mutation_count_sigprofiler,
                            treatment_mutation_count_siganalyzer
                        ])
                        mean_exposure_methods = np.mean([
                            proportion_treatment_mutations_sigprofiler,
                            proportion_treatment_mutations_siganalyzer
                        ])

                        coherent_samples_ttype[drug][dic_primary_full[sample]]['exposure'].append(
                            mean_exposure_methods
                        )
                        coherent_samples_ttype[drug][dic_primary_full[sample]]['count'].append(
                            mean_count_methods
                        )

                        coherent_samples_samplelevel[drug][sample] = mean_count_methods

                    # create a list because this sample will be forbidden in some of the analyses (only aging)
                    elif (class_sample == 'not-concordant') & (drug != 'AGING'):
                        forbidden_samples.add(sample)

    wanted_ttypes = ['Colon-Rectum', 'Esophagus', 'Lung', 'Urinary-tract', 'Uterus', 'Ovary', 'Breast']
    # Do Sigmoid plots
    (samples_mutations_treatment_coherent, samples_mutations_treatment_sigpro,
     samples_mutations_treatment_sigan) = sigmoid_plots(
        wanted_ttypes,
        allowed_samples,
        forbidden_samples,
        dic_primary_full,
        df_siganalyzer,
        df_sigprofiler,
        siganalyzer_sigs,
        sigprofiler_sigs,
        total_samples,
        list_samples_sigprofiler,
        treated_samples_dictionary
    )
    # do wiskeys plots
    plot_wiskeys(wanted_ttypes, coherent_samples_ttype, 'coherent')
    plot_wiskeys(wanted_ttypes, exposures_drug_siganalyzer_ttype, 'siganalyzer')
    plot_wiskeys(wanted_ttypes, exposures_drug_sigprofiler_ttype, 'sigprofiler')

    # do barplots from the main figure
    plot_barplots(wanted_ttypes, class_samples_according_coherence)

    write_supptable(
        df_siganalyzer, wanted_ttypes, forbidden_samples, coherent_samples_samplelevel,
        exposures_drug_sigprofiler, exposures_drug_siganalyzer, samples_mutations_treatment_coherent,
        samples_mutations_treatment_sigpro, samples_mutations_treatment_sigan, dic_primary_full
    )

    return coherent_samples_samplelevel, exposures_drug_siganalyzer, exposures_drug_sigprofiler


# this will create all aging-related figures
def aging_plots(coherent_samples_samplelevel, exposures_drug_siganalyzer, exposures_drug_sigprofiler, dic_primary_full):

    counter_length = pickle.load(gzip.open('data/clinical_data/dic_counter_total_days_treatment.pckl.gz'))
    prob_ttypes = pickle.load(gzip.open("data/probability_MB/results_probability_sigpro_20190502.pckl.gz"))

    dic_age = return_aging()
    dic_results, dic_results_ttype = proportion_aging_signature(
        coherent_samples_samplelevel, counter_length, dic_age, dic_primary_full
    )

    # FIG6a
    plot_aging_twosigs(dic_results, "coherent")
    # Fig6B
    plot_lines_probabilities(prob_ttypes, "coherent")

    # Fig6c and supps
    plot_lines_time(prob_ttypes, dic_results_ttype, "gene", "coherent")
    plot_lines_time(prob_ttypes, dic_results_ttype, "cgc", "coherent")

    ########
    # EDF10
    ########

    # --------
    # Restricted list
    # --------
    if os.path.isfile('data/clinical_data/sample_ID_martijn.txt'):
        samples = pd.read_csv('data/clinical_data/sample_ID_martijn.txt', sep='\t', names=['sample'])
        samples_list = samples['sample'].tolist()
        counter_length_martijn = {k: counter_length[k] for k in counter_length if k in samples_list}

        new_exposures = defaultdict(lambda: defaultdict(int))
        for treatment, dt in coherent_samples_samplelevel.items():
            for sample, v in dt.items():
                if sample in samples_list:
                    new_exposures[treatment][sample] = v

        dic_results, dic_results_ttype = proportion_aging_signature(
            new_exposures, counter_length_martijn, dic_age, dic_primary_full
        )
        plot_aging_twosigs(dic_results, "restricted")
        plot_lines_time(prob_ttypes, dic_results_ttype, "gene", "restricted")
        plot_lines_time(prob_ttypes, dic_results_ttype, "cgc", "restricted")

    # ---------
    # literature
    # ---------
    counter_literature = get_complete_regime()
    dic_results, dic_results_ttype = proportion_aging_signature(
        coherent_samples_samplelevel, counter_literature, dic_age, dic_primary_full
    )

    plot_aging_twosigs(dic_results, "literature")
    plot_lines_probabilities(prob_ttypes, "literature")
    plot_lines_time(prob_ttypes, dic_results_ttype, "gene", "literature")
    plot_lines_time(prob_ttypes, dic_results_ttype, "cgc", "literature")

    dic_results, dic_results_ttype = proportion_aging_signature(
        exposures_drug_siganalyzer, counter_length, dic_age, dic_primary_full
    )
    # EDF8d
    plot_aging_twosigs(dic_results, "SignatureAnalyzer")
    dic_results, dic_results_ttype = proportion_aging_signature(
        exposures_drug_sigprofiler, counter_length, dic_age, dic_primary_full
    )
    # EDF9d
    plot_aging_twosigs(dic_results, "SigProfiler")


def generate_figures_coherent_methods():

    path_base = 'data/hartwig/signatures/extraction/results/{}/snvs/exposures/'
    path_siganalyzer = path_base.format('SignatureAnalyzer') + "Pan_full/Pan_full.exposures.tsv"
    path_sigprofiler = path_base.format('SigProfiler') + "PanNoSkinNoUnknown.snvs/PanNoSkinNoUnknown.snvs.exposures.tsv"

    df_siganalyzer = read_exposures(path_siganalyzer)
    df_sigprofiler = read_exposures(path_sigprofiler)

    treated_samples_dictionary = pickle.load(gzip.open('data/clinical_data/hartwig_FDA_treatments_specific.pckl.gz'))

    dic_primary_full, _ = return_metadata()

    siganalyzer_sigs = {
        'CISPLATIN': ['21_SBS31_0.953955_1', '14_1'],
        'CARBOPLATIN': ['21_SBS31_0.953955_1', '25_1'],
        '5-FU_CAPE': ['31_SBS17b_0.968799_1'],
        'OXALIPLATIN': ['14_1', '37_1'],
        'AGING': ['17_SBS1_0.995019_1']
    }
    sigprofiler_sigs = {
        'CISPLATIN': ['1_SBS31_0.968153_0.98'],
        'CARBOPLATIN': ['1_SBS31_0.968153_0.98'],
        'OXALIPLATIN': ['20_0.92'],
        '5-FU_CAPE': ['19_SBS17b_0.961548_0.99'],
        'AGING': ['23_SBS1_0.996494_1.0'],
    }

    coherent_samples_samplelevel, exposures_drug_siganalyzer, exposures_drug_sigprofiler = analysis_coherence(df_siganalyzer, df_sigprofiler, siganalyzer_sigs, sigprofiler_sigs, treated_samples_dictionary, dic_primary_full)
    aging_plots(coherent_samples_samplelevel, exposures_drug_siganalyzer, exposures_drug_sigprofiler, dic_primary_full)


if __name__ == '__main__':
    generate_figures_coherent_methods()
