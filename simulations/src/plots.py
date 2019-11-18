""" Plotting functions to generate figures panels """

import os
import json
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import seaborn

from plot_signatures import plot_signature
from recovery import signature_recovery, exposure_recovery, exposure_recovery_nontreated, concordance_cc


def params_plot():

    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'


params_plot()


# jitter for plots

def jitter(arr):

    std = .01 * (max(arr) - min(arr))
    return arr + np.random.randn(len(arr)) * std


# fit in boxplots

def boxplot_fit(x_arr, y_arr):

    x_values = sorted(list(set(x_arr)))
    y_boxed = [[] for _ in x_values]
    for i, y in enumerate(y_arr):
        j = x_values.index(x_arr[i])
        y_boxed[j].append(y)
    return y_boxed


# plot profile

def plot_profiles(sig, signatures_dict, **kwargs):

    fig, ax = plt.subplots(figsize=(20, 5))

    plot_signature(signatures_dict[sig], ax=ax, **kwargs)
    ax.set_title(sig)
    fig.savefig(f'./figures/profile_{sig}.png', dpi=200, bbox_inches='tight')


# plot burden distribution

def plot_burden_distribution(breast, breastlungcolon):

    fig, ax = plt.subplots(figsize=(15, 5), ncols=2)

    breast_burden = np.log10(breast.sum(axis=0).values)
    breastlungcolon_burden = np.log10(breastlungcolon.sum(axis=0).values)
    seaborn.distplot(breast_burden, ax=ax[0])
    seaborn.distplot(breastlungcolon_burden, ax=ax[1])

    # axis 0

    ax[0].set_title('Breast')
    xticks = [2, 3, 4, 5, 6, 7]
    ax[0].set_xticks(xticks)
    ax[0].set_xticklabels(labels=ax[0].set_xticklabels(labels=list(map(lambda x: f'$10^{str(x)}$', xticks))))
    ax[0].set_xlabel('no. mutations')
    ax[0].set_ylabel('density')
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)

    # axis 1

    ax[1].set_xticks(xticks)
    ax[1].set_xticklabels(labels=ax[1].set_xticklabels(labels=list(map(lambda x: f'$10^{str(x)}$', xticks))))
    ax[1].set_title('Breast-Lung-Colorectal')
    ax[1].set_xlabel('no. mutations')
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)

    plt.savefig('./figures/burden_distribution.png', dpi=200, bbox_inches='tight')

    plt.show()


# plot all cosine similarities

def sig_scatter_cosine(sig, signatures, offset, color, folder, ax):

    cosines = []
    x = []
    twin_signals = []
    major_labels = [10, 25, 50, 100, 150]  # representing n_treated
    for i, n_treated in enumerate(major_labels):
        cos = list(map(lambda x: signature_recovery(sig, n_treated, x, signatures, folder), range(25)))
        cosines += list(zip(*cos))[1]
        twin_signals += list(zip(*cos))[2]
        x += [i + offset] * len(cos)
    colors = ['red' if s else color for s in twin_signals]

    # boxplot layer

    alpha_param = 0.3
    y_boxplot = boxplot_fit(x, cosines)
    ax.boxplot(y_boxplot, positions=np.arange(5)+offset, showfliers=False, boxprops={'alpha': alpha_param},
               whiskerprops={'alpha': alpha_param}, capprops={'alpha': alpha_param}, widths=0.3,
               medianprops={'color': 'black', 'alpha': alpha_param})

    # scatter plot layer

    ax.scatter(jitter(x), cosines, color=colors, alpha=0.4, s=10, label=sig)


def scatter_cosine(sigs, signatures, folder, title):

    fig, ax = plt.subplots(figsize=(8, 5))

    sig_scatter_cosine(sigs[0], signatures, -0.25, 'blue', folder, ax)
    sig_scatter_cosine(sigs[1], signatures, 0.1, 'maroon', folder, ax)

    # plot layout

    major_labels = [10, 25, 50, 100, 150]  # representing n_treated
    major_ticks = np.arange(5)
    ax.set_ylim(0.67, 1.03)
    ax.set_xticks(ticks=major_ticks)
    ax.set_xticklabels(labels=map(str, major_labels))
    ax.set_xlabel('no. samples exposed')
    ax.set_ylabel('cosine similarity')
    ax.set_title(f'Similarity of Reconstructed vs. Injected ({title})')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(loc=4)

    cohort = os.path.basename(folder)
    fig.savefig(f'./figures/{cohort}.combined.cosine.svg', dpi=400, bbox_inches='tight')
    plt.show()


"""
plot signature-treatment dependencies: 
- log-odds-ratio 
- significance 
- effect size (log-fold-change)
"""


def plot_lor(sig, n_treated, replicate, signatures, regression_folder, deconstruction_folder, shuffle=None):

    index, cosine, twin_signal = signature_recovery(sig, n_treated, replicate, signatures, deconstruction_folder)

    if shuffle is not None:
        fn = os.path.join(regression_folder, f'shuffle.{shuffle}.regression_lor.{sig}.{n_treated}.{replicate}.json')
    else:
        fn = os.path.join(regression_folder, f'regression_lor.{sig}.{n_treated}.{replicate}.json')

    with open(fn) as f:
        merge_res = json.load(f)

    pval_red, pval_pink, pval_other = [], [], []
    for k in sorted(merge_res, key=lambda v: np.mean(merge_res[v]), reverse=True):
        x, y = [], []
        for lor in merge_res[k]:
            x.append(k)
            y.append(lor)
        pval = np.sum(np.array(y) < 0) / len(y)
        if k == str(index[0]):
            pval_red.append(pval)
        elif k in list(map(str, index)):
            pval_pink.append(pval)
        else:
            pval_other.append(pval)

    return pval_red, pval_pink, pval_other


def plot_all_lor(sig, n_treated, signatures, regression_folder, deconstruction_folder, shuffle=None):

    pvals_red, pvals_pink, pvals_other = [], [], []
    for replicate in range(25):
        p_red, p_pink, p_other = plot_lor(sig, n_treated, replicate, signatures,
                                          regression_folder, deconstruction_folder,
                                          shuffle=shuffle)
        pvals_red += p_red
        pvals_pink += p_pink
        pvals_other += p_other

    return pvals_red, pvals_pink, pvals_other


def lor_mean(sig, n_treated, replicate, signatures, regression_folder, deconstruction_folder, shuffle=None):

    index, cosine, twin_signal = signature_recovery(sig, n_treated, replicate, signatures, deconstruction_folder)

    if shuffle is not None:
        fn = os.path.join(regression_folder,
                          f'shuffle.{shuffle}.regression_lor.{sig}.{n_treated}.{replicate}.json')
    else:
        fn = os.path.join(regression_folder,
                          f'regression_lor.{sig}.{n_treated}.{replicate}.json')

    with open(fn, 'rt') as f:
        merge_res = json.load(f)

    mean_lor, x = [], []
    for k in sorted(merge_res, key=lambda v: np.mean(merge_res[v]), reverse=True):
        if k in list(map(str, index)):
            mean_lor.append(np.mean(merge_res[k]))
            x.append(0)
        else:
            mean_lor.append(np.mean(merge_res[k]))
            x.append(1)

    return x, mean_lor


def plot_lor_mean(sig, signatures, regression_folder, deconstruction_folder, shuffle=None):

    fig, ax = plt.subplots(figsize=(8, 5))

    all_x, values = [], []
    for i, n_treated in enumerate([10, 25, 50, 100, 150]):
        for replicate in range(25):
            try:
                x, mean_lor = lor_mean(sig, n_treated, replicate, signatures,
                                       regression_folder, deconstruction_folder,
                                       shuffle=shuffle)
                all_x += list((np.array(x) - 0.5) * (0.2 / 0.5) + i)
                values += list(mean_lor)
            except FileNotFoundError:
                continue

    ax.scatter(jitter(all_x), values, alpha=0.4, color='black', s=10)
    ax.hlines(0, -0.4, 4.4, linestyles='dashed', color='red', alpha=0.5)
    ax.set_ylabel('mean log odds ratio')
    ax.tick_params(length=0, which='major', pad=65, labelsize=14, axis='x')
    ax.tick_params(length=0, which='minor', pad=5, labelsize=14, rotation=90)
    ax.set_xticks(sum([[i - 0.2, i + 0.2] for i in range(5)], []), minor=True)
    ax.set_xticklabels([sig, 'other'] * 5, minor=True)
    ax.set_xticks(range(5), minor=False)
    ax.set_xticklabels(list(map(lambda x: f'{x} exposed', [10, 25, 50, 100, 150])))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    cohort = os.path.basename(deconstruction_folder)
    if shuffle is not None:
        fn = f'figures/{cohort}.shuffle.{shuffle}.{sig}.mean_logodds_plot.png'
        title = f'{cohort} : {sig} : {shuffle}% shuffled'
    else:
        fn = f'figures/{cohort}.{sig}.mean_logodds_plot.png'
        title = f'{cohort} : {sig}'

    ax.set_title(title)
    fig.savefig(fn, dpi=200, bbox_inches='tight')
    plt.show()


def plot_lor_significance(sig, signatures, regression_folder, deconstruction_folder, shuffle=None):

    fig, ax = plt.subplots(figsize=(8, 5))

    all_x, values, colors = [], [], []
    for i, n_treated in enumerate([10, 25, 50, 100, 150]):
        x = []
        try:
            pvals_red, pvals_pink, pvals_other = plot_all_lor(sig, n_treated, signatures,
                                                              regression_folder, deconstruction_folder,
                                                              shuffle=shuffle)
        except FileNotFoundError:
            continue

        values += list(map(lambda x: -np.log10(x + 1e-4), pvals_red))
        x += [0] * len(pvals_red)
        colors += ['red'] * len(pvals_red)
        values += list(map(lambda x: -np.log10(x + 1e-4), pvals_pink + pvals_other))
        x += [1] * len(pvals_pink)
        colors += ['pink'] * len(pvals_pink)
        x += [1] * len(pvals_other)
        colors += ['black'] * len(pvals_other)
        all_x += list((np.array(x) - 0.5) * (0.2 / 0.5) + i)

    ax.scatter(jitter(all_x), values, alpha=0.4, color=colors, s=10)
    ax.hlines(3, -0.4, 4.4, linestyles='dashed', color='red', alpha=0.5)
    ax.set_ylabel('-log pvalue')
    ax.tick_params(length=0, which='major', pad=50, labelsize=14, axis='x')
    ax.tick_params(length=0, which='minor', pad=5, labelsize=14, rotation=90)
    ax.set_xticks(sum([[i - 0.2, i + 0.2] for i in range(5)], []), minor=True)
    ax.set_xticklabels(['true', 'false'] * 5, minor=True)
    ax.set_xticks(range(5), minor=False)
    ax.set_xticklabels(list(map(lambda x: f'{x} exposed', [10, 25, 50, 100, 150])))
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_yticklabels([0, 1, 2, 3, '$\infty$'])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    cohort = os.path.basename(deconstruction_folder)
    if shuffle is not None:
        fn = f'figures/{cohort}.shuffle.{shuffle}.{sig}.significance_plot.png'
        title = f'{cohort} : {sig} : {shuffle}% shuffled'
    else:
        fn = f'figures/{cohort}.{sig}.significance_plot.png'
        title = f'{cohort} : {sig}'

    ax.set_title(title)
    fig.savefig(fn, dpi=200, bbox_inches='tight')
    plt.show()


def fold_change(sig, n_treated, replicate, signatures, regression_folder, deconstruction_folder, shuffle=None):

    index, cosine, twin_signal = signature_recovery(sig, n_treated, replicate, signatures, deconstruction_folder)

    if shuffle is not None:
        fn = os.path.join(regression_folder,
                          f'shuffle.{shuffle}.regression_summary.{sig}.{n_treated}.{replicate}.tsv')
    else:
        fn = os.path.join(regression_folder,
                          f'regression_summary.{sig}.{n_treated}.{replicate}.tsv')

    df = pd.read_csv(fn, sep='\t')

    x, y = [], []
    x.append(0)
    y.append(df.loc[index[0], 'effect_size'])
    x += [1 for i in df.loc[:, 'signature'] if i not in [index[0]]]
    y += list(df[~df.signature.isin([index[0]])]['effect_size'].values)

    return x, y


def plot_fold_change(sig, signatures, regression_folder, deconstruction_folder, shuffle=None, violin=False):

    fig, ax = plt.subplots(figsize=(10, 3))

    ax.tick_params(length=0, which='major', pad=60, labelsize=14, axis='x')
    ax.tick_params(length=0, which='minor', pad=5, labelsize=14, rotation=90)
    ax.set_xticks(sum([[i - 0.2, i + 0.2] for i in range(5)], []), minor=True)
    ax.set_xticklabels([sig, 'other'] * 5, minor=True)
    ax.set_xticks(range(5), minor=False)
    ax.set_xticklabels(list(map(lambda x: f'treated={x}', [10, 25, 50, 100, 150])))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    for i, n_treated in enumerate([10, 25, 50, 100, 150]):
        for replicate in range(1, 11):
            try:
                x, y = fold_change(sig, n_treated, replicate, signatures,
                                   regression_folder, deconstruction_folder,
                                   shuffle=shuffle)
                if violin:
                    seaborn.violinplot((np.array(x) - 0.5) * (0.2 / 0.5) + i, y,
                                       alpha=0.3, color='#DCDCDC', scale='count', inner=None, ax=ax)
                ax.scatter(jitter((np.array(x) - 0.5) * (0.2 / 0.5) + i), y,
                           alpha=0.4, color='black', s=10)
            except FileNotFoundError:
                continue

    ax.hlines(2, -0.4, 4.4, linestyles='dashed', color='red', alpha=0.4)
    ax.set_ylabel('fold change')

    cohort = os.path.basename(deconstruction_folder)
    if shuffle is not None:
        title = f'{cohort} : {sig} : {shuffle}% shuffle'
        fn = f'figures/{cohort}.shuffle.{shuffle}.{sig}.fold_change_plot.png'
    else:
        title = f'{cohort} : {sig}'
        fn = f'figures/{cohort}.{sig}.fold_change_plot.png'
    ax.set_title(title)
    fig.savefig(fn, dpi=200, bbox_inches='tight')
    plt.show()


# volcano plots

def get_effects_pvals(sig, n_treated, replicate, signatures, regression_folder, deconstruction_folder, shuffle=None):

    index, cosine, twin_signal = signature_recovery(sig, n_treated, replicate, signatures, deconstruction_folder)
    if shuffle is not None:
        fn = os.path.join(regression_folder,
                          f'shuffle.{shuffle}.regression_lor.{sig}.{n_treated}.{replicate}.json')
        with open(fn, 'rt') as f:
            merge_res = json.load(f)
        summary = pd.read_csv(os.path.join(regression_folder,
                                           f'shuffle.{shuffle}.regression_summary.{sig}.{n_treated}.{replicate}.tsv'),
                              sep='\t')
    else:
        fn = os.path.join(regression_folder, f'regression_lor.{sig}.{n_treated}.{replicate}.json')
        with open(fn, 'rt') as f:
            merge_res = json.load(f)
        summary = pd.read_csv(os.path.join(regression_folder,
                                           f'regression_summary.{sig}.{n_treated}.{replicate}.tsv'),
                              sep='\t')

    effects, pvals, colors = [], [], []
    for k in sorted(merge_res, key=lambda v: np.mean(merge_res[v]), reverse=True):
        y = []
        for lor in merge_res[k]:
            y.append(lor)
        pval = np.sum(np.array(y) < 0) / len(y)
        pvals.append(pval)
        if int(k) == index[0]:
            colors.append('red')
        else:
            colors.append('black')
        effects.append(summary.loc[int(k), 'effect_size'])

    return effects, pvals, colors


def plot_volcano(sig, n_treated, signatures, regression_folder, deconstruction_folder, ax, shuffle=None):

    x, y, c = [], [], []
    for replicate in range(25):
        try:
            effects, pvals, colors = get_effects_pvals(sig, n_treated, replicate, signatures,
                                                       regression_folder, deconstruction_folder,
                                                       shuffle=shuffle)
            x += effects
            y += list(map(lambda x: -np.log10(x + 1e-4), pvals))
            c += colors
        except FileNotFoundError:
            continue

    ax.scatter(x, y, alpha=0.4, color=c, s=10)
    ax.hlines(3, -0.1, np.max(x), linestyles='dashed', color='red', alpha=0.5)
    ax.vlines(2, -0.1, np.max(y), linestyles='dashed', color='red', alpha=0.5)
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_yticklabels([0, 1, 2, 3, '$\infty$'])
    ax.set_ylabel('$-\log_{10}$ p-value')
    ax.set_xlabel('effect size (fold change)')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    cohort = os.path.basename(deconstruction_folder)
    if shuffle is not None:
        title = f'{cohort} : {sig}\n{n_treated} exposed : {shuffle}% shuffled'
    else:
        title = f'{cohort} : {sig}\n{n_treated} exposed'
    ax.set_title(title)


def plot_volcano_panel(sig, signatures, regression_folder, deconstruction_folder, shuffle=None):

    fig, ax = plt.subplots(figsize=(25, 5), ncols=5)

    for i, n_treated in enumerate([10, 25, 50, 100, 150]):
        plot_volcano(sig, n_treated, signatures, regression_folder, deconstruction_folder, ax[i], shuffle=shuffle)

    cohort = os.path.basename(deconstruction_folder)
    if shuffle is not None:
        fn = f'figures/{cohort}.shuffle.{shuffle}.{sig}.volcano.png'
    else:
        fn = f'figures/{cohort}.{sig}.volcano.png'

    fig.savefig(fn, dpi=200, bbox_inches='tight')
    plt.show()


# plot results of cotreatment analysis

def plot_cotreatment(sig, cotreatment_folder, symmetric, one_to_one):

    fig, ax = plt.subplots()

    all_x, all_y, overlaps = [], [], []
    symmetric_label = 'symmetric' if symmetric else 'asymmetric'
    one_to_one_label = 'one_to_one' if one_to_one else 'many_to_one'

    for overlap in [25, 50, 75]:
        try:
            fn = os.path.join(cotreatment_folder,
                              f'overlap.{overlap}.{sig}.{symmetric_label}.{one_to_one_label}'
                              f'.cotreatment_fold_change.json')
            with open(fn, 'rt') as f:
                fold_change_value = json.load(f)
            x, y = list(zip(*[(i, v) for i, fc in enumerate(fold_change_value) for v in fc]))
            all_x += x
            all_y += y
            overlaps += [overlap] * len(x)
        except FileNotFoundError:
            continue

    data = pd.DataFrame({'x': all_x, 'y': all_y, 'overlap': overlaps})
    seaborn.swarmplot('x', 'y', hue='overlap', data=data, dodge=True, s=3, alpha=1)
    ax.set_xticklabels(labels=list(map(str, [10, 25, 50, 100, 150])))
    ax.set_xlabel('no. samples exposed')
    ax.set_ylabel('$\log_{10}$ fold change')
    ax.set_ylim(-5, 5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.hlines(0, -0.5, 4.5, color='red', linestyles='dashed')
    ax.legend(title='Overlap', loc=4, labels=['25%', '50%', '75%'])

    cohort = os.path.basename(cotreatment_folder)
    ax.set_title(f'{cohort} : {sig} : {symmetric_label} : {one_to_one_label}')
    plt.savefig(f'figures/{cohort}.{sig}.{symmetric_label}.{one_to_one_label}.cotreatments.png',
                dpi=200, bbox_inches='tight')


# exposure recovery plot

def plot_exposure_recovery(sig, tumor_type, n_treated, signatures, catalogue_folder, deconstruction_folder):

    fig, ax = plt.subplots(figsize=(8, 5))

    for replicate in range(25):
        try:
            burden, exposure = exposure_recovery(sig, n_treated, replicate, signatures,
                                                 catalogue_folder, deconstruction_folder)
        except FileNotFoundError:
            continue

        ax.scatter(burden, exposure, alpha=0.4, color='black', s=15)

    ax.set_xlabel('synthetic')
    ax.set_ylabel('reconstructed')
    ax.set_ylim(0, 5000)
    ax.set_xlim(0, 5000)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    x = np.linspace(0, 5000, 5000)
    ax.plot(x, x, '--', color='red', alpha=0.4)

    cohort = os.path.basename(deconstruction_folder)
    ax.set_title(f'Exposure Recovery\n{tumor_type} : {sig} : {n_treated} samples exposed')
    fig.savefig(f'figures/{cohort}.exposure_recovery.{sig}.{n_treated}.png', dpi=200, bbox_inches='tight')


def scatter_error(sig, signatures, catalogue_folder, deconstruction_folder, option='error', hline=0, violin=False):
    """
    option:
        'error',
        'correlation',
        'relative_residual',
        'concordance'
    """

    fig, ax = plt.subplots()

    # plot layout

    major_labels = [10, 25, 50, 100, 150]
    if option == 'error':
        ax.set_ylabel('mean relative error (%)')
    elif option == 'correlation':
        ax.set_ylabel('correlation')
    elif option == 'relative_residual':
        ax.set_ylabel('relative residual')
    elif option == 'concordance':
        ax.set_ylabel('concordance')

    # scatter plot

    values = []
    x = []
    for i, n_treated in enumerate(major_labels):
        for replicate in range(25):
            burden, exposure = exposure_recovery(sig, n_treated, replicate, signatures,
                                                 catalogue_folder, deconstruction_folder)
            if option == 'error':
                value = np.mean(np.abs((burden - exposure) / burden))
            elif option == 'correlation':
                value = [scipy.stats.pearsonr(burden, exposure)[0]]
                hline = 1
            elif option == 'relative_residual':
                value = (burden - exposure) / burden
            elif option == 'concordance':
                value = [concordance_cc(exposure, burden)]
                hline = 1
            x += [i] * len(value)
            values += list(value)
    if violin:
        seaborn.violinplot(x, values, alpha=0.3, color='#DCDCDC', scale='count', inner=None, ax=ax)
    ax.scatter(jitter(x), values, alpha=0.1, color='black', s=10)

    major_ticks = np.arange(5)
    if option == 'error':
        ax.set_yticks([0, 1, 2, 3])
        ax.set_yticklabels(labels=map(str, [0, 100, 200, 300]))
    ax.set_xticks(ticks=major_ticks)
    ax.set_xticklabels(labels=map(str, major_labels))
    ax.set_title(sig)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.hlines(hline, -0.5, 4.5, color='red', linestyles='dashed')
    ax.set_xlabel('no. exposed samples')

    cohort = os.path.basename(deconstruction_folder)
    ax.set_title(f'{cohort} : {sig}')
    fig.savefig(f'figures/{cohort}.{sig}.exposure_{option}.png', dpi=200, bbox_inches='tight')
    plt.show()


def average_exposure_vs_average_injected_nontreated(sig, tumor_type, signatures, catalogue_folder, deconstruction_folder, violin=False):

    fig, ax = plt.subplots()

    x = []
    values = []
    major_labels = [10, 25, 50, 100, 150]
    for i, n_treated in enumerate(major_labels):
        for replicate in range(1, 11):

            exposure = exposure_recovery_nontreated(sig, n_treated, replicate, signatures,
                                                    catalogue_folder, deconstruction_folder)
            try:
                # load synthetic burden

                burden_fn = f'{sig}.{n_treated}.{replicate}.burden.json'
                with open(os.path.join(catalogue_folder, burden_fn), 'rt') as f:
                    burden = json.load(f)
                burden = burden[:n_treated]
                x += [i]
                values.append(np.mean(exposure) / np.mean(burden))
            except FileNotFoundError:
                continue

    if violin:
        seaborn.violinplot(x, values, alpha=0.3, color='#DCDCDC', scale='count', inner=None, ax=ax)
    ax.scatter(jitter(x), values, alpha=0.4, color='black', s=10)

    major_ticks = np.arange(5)
    ax.set_xticks(ticks=major_ticks)
    ax.set_xticklabels(labels=map(str, major_labels))
    ax.set_xlabel('no. exposed samples')
    ax.set_ylabel('average exposure over\naverage injected burden')
    ax.set_ylim(-0.08, 0.6)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.hlines(0, -0.5, 4.5, color='red', linestyles='dashed')

    cohort = os.path.basename(deconstruction_folder)
    ax.set_title(f'{tumor_type} : {sig}')
    fig.savefig(f'figures/{cohort}.{sig}.exposure_recovery_notreated.png', dpi=200, bbox_inches='tight')
