import os
import sys
from collections import Counter

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from statsmodels.stats import proportion

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
scripts_path = os.path.abspath(os.path.join(__file__, '..', '..',))
sys.path.insert(0, scripts_path)

from utils.orders import order_muts, order_to_plot_dbs, order_to_plot_indel, order_to_plot_snvs
from utils.plots_utils import config_params


def poisson_exact(y1, y2, r_d=1, alternative='two-sided'):

    bp = r_d / (1 + r_d)
    y_total = y1 + y2
    stat = None
    pvalue = proportion.binom_test(y1, y_total, prop=bp, alternative=alternative)

    return pvalue


# split into even chunks
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


# function to plot the SNV processes
def plot_snvs(sig, title, outpath, ttype):

    config_params(3)

    fig, axs = plt.subplots(
        nrows=2, ncols=1, figsize=(3.2, 1), gridspec_kw={'height_ratios': [1, 9]}
    )
    order_plot = order_muts("snv")

    vals = []
    colors = []
    colors_mut = [
        '#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5'
    ]
    bot = -0.5
    for ix, c in enumerate(chunks(sig, 16)):
        colors.extend([colors_mut[ix] for s in c])
        axs[0].barh(1, 16, left=bot, color=colors_mut[ix])
        bot += 16
        vals.extend(c)

    axs[0].set_xlim(-1, 96)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    x = [i for i in range(len(vals))]

    axs[1].axhline(y=0.05, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.1, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.15, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)

    axs[1].bar(x, vals, color=colors, width=0.8, linewidth=0, align='center')
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(
        ['{}{}{}'.format(a[0], a[2], a[-1]) for a in order_plot],
        verticalalignment="center", ha='center', rotation=90, fontsize=2,
        color='grey'
    )

    plt.tight_layout()
    plt.xlim(-1, 96)

    axs[1].spines['top'].set_visible(False)
    axs[1].set_ylabel('Relative Probability')

    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=0.5)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.tick_params(axis='both', which='both', bottom=False, left=False)
    plt.savefig('{}/processes/{}/{}.{}.png'.format(outpath, ttype, ttype, title), dpi=300, bbox_inches='tight')
    plt.savefig('{}/processes/{}/{}.{}.svg'.format(outpath, ttype, ttype, title))

    plt.close()


# function to plot the DBS processes
def plot_dbs(sig, title, outpath, ttype):
    config_params(3)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(3.2, 1),
                            gridspec_kw={'height_ratios': [1, 9]})

    vals = []
    colors_mut = ['#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5']

    dbs_color = {
        'AC': '#a6cee3', 'AT': '#1f78b4', 'CC': '#b2df8a', 'CG': '#33a02c', 'CT': '#fb9a99',
        'GC': '#e3211d', 'TA': '#fdbf6f', 'TC': '#ff7f00', 'TG': '#cab2d6',
        'TT': '#6a3d9a',
    }

    order_color = [
        '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
        '#e3211d', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a'
    ]

    order_dbs_list = order_to_plot_dbs()

    vals = sig

    colors = [dbs_color[db.split('_')[0]] for db in order_dbs_list]
    counter_colors = Counter(colors)

    bot = -0.5

    for c in order_color:
        axs[0].barh(1, counter_colors[c], left=bot, color=c, align='center',)
        bot += counter_colors[c]

    axs[0].spines['top'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    axs[0].set_xlim(-1, 78)

    x = [i for i in range(len(vals))]

    axs[1].axhline(y=0.05, xmin=-1, xmax=96, lw=0.3, color='grey', alpha=0.2)
    axs[1].axhline(y=0.1, xmin=-1, xmax=96, lw=0.3, color='grey', alpha=0.2)
    axs[1].axhline(y=0.15, xmin=-1, xmax=96, lw=0.3, color='grey', alpha=0.2)

    axs[1].bar(x, vals, color=colors, width=0.8, linewidth=0, align='center', alpha=1)
    axs[1].set_xticks(x)

    axs[1].set_xticklabels(['{}{}'.format(a[-2], a[-1],) for a in order_dbs_list],
                           rotation=90, fontsize=2, verticalalignment="center", ha='center',
                           color='grey')

    axs[1].set_xlim(-1, 78)

    axs[1].spines['top'].set_visible(False)
    axs[1].set_ylabel('Relative Probability')

    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=0.5)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.tick_params(axis='both', which='both', bottom=False, left=False)
    plt.tight_layout()
    plt.savefig('{}/processes/{}/{}.{}.png'.format(outpath, ttype, ttype, title), dpi=300, bbox_inches='tight')
    plt.savefig('{}/processes/{}/{}.{}.svg'.format(outpath, ttype, ttype, title))

    plt.close()


# function to plot the indel processes
def plot_indel(sig, title, outpath, ttype):

    config_params(3)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(3.2, 1),
                            gridspec_kw={'height_ratios': [1, 9]})

    vals = sig

    colors = ['#fdbe6f'] * 6 + ['#ff8002'] * 6 + ['#b0dd8b'] * 6 + ['#36a12e'] * 6 + \
             ['#fdcab5'] * 6 + ['#fc8a6a'] * 6 + ['#f14432'] * 6 + ['#bc191a'] * 6 + \
             ['#d0e1f2'] * 6 + ['#94c4df'] * 6 + ['#4a98c9'] * 6 + ['#1764ab'] * 6 + \
             ['#e1e1ef'] * 1 + ['#b6b6d8'] * 2 + ['#8683bd'] * 3 + ['#62409b'] * 5

    order_colors = ['#fdbe6f', '#ff8002', '#b0dd8b', '#36a12e', '#fdcab5', '#fc8a6a', '#f14432',
                    '#bc191a', '#d0e1f2', '#94c4df', '#4a98c9', '#1764ab', '#e1e1ef', '#b6b6d8', '#8683bd',
                    '#62409b']

    counter_colors = Counter(colors)

    bot = -0.5

    for c in order_colors:
        axs[0].barh(1, counter_colors[c], left=bot, color=c)
        bot += counter_colors[c]

    axs[0].spines['top'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    axs[0].set_xlim(-1, 83)

    x = [i for i in range(len(vals))]

    axs[1].axhline(y=0.05, xmin=-1, xmax=83, lw=0.3, color='grey', alpha=0.3)
    axs[1].axhline(y=0.1, xmin=-1, xmax=83, lw=0.3, color='grey', alpha=0.3)
    axs[1].axhline(y=0.15, xmin=-1, xmax=83, lw=0.3,  color='grey', alpha=0.3)
    axs[1].bar(x, vals, color=colors, width=0.7, linewidth=0, align='center', alpha=1)
    axs[1].set_xticks(x)
    axs[1].set_xticklabels([i.split('_')[-1] for i in order_to_plot_indel()], fontsize=2,
                           verticalalignment="center", ha='center', color='grey')
    axs[1].set_xlim(-1, 83)

    plt.tight_layout()

    axs[1].spines['top'].set_visible(False)
    axs[1].set_ylabel('Relative Probability')

    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=0.5)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    plt.tick_params(axis='both', which='both', bottom=False, left=False)
    plt.savefig('{}/processes/{}/{}.{}.png'.format(outpath, ttype, ttype, title), dpi=300, bbox_inches='tight')
    plt.savefig('{}/processes/{}/{}.{}.svg'.format(outpath, ttype, ttype, title))

    plt.close()


# function to plot snvs in bias
def plot_bias_snv(d_clustered, d_unclustered, outpath, ttype, label1, label2):

    config_params(3)

    final_order = []
    order_plot = order_to_plot_snvs()

    for o in order_plot:
        final_order.append('{}_{}'.format(o, label1))
        final_order.append('{}_{}'.format(o, label2))

    for s in d_clustered.columns:

        sig = d_clustered[s]

        fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10, 4))

        colors = []
        colors_mut = ['#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5']

        for ix, c in enumerate(chunks(sig, 16)):
            colors.extend([colors_mut[ix] for s in c])
            colors.extend([colors_mut[ix] for s in c])

        sig = []
        for lag, lead in zip(d_clustered[s], d_unclustered[s]):
            sig.append(lag)
            sig.append(lead)

        # get significance per type
        toplot_full = []
        colors_joinbar = []
        marks = []
        for ix, c in enumerate(chunks(sig, 32)):

            c_condition1 = np.sum(c[::2])
            c_condition2 = np.sum(c[1::2])
            colors_joinbar.append(colors_mut[ix])
            colors_joinbar.append(colors_mut[ix])
            toplot_full.append(c_condition1)
            toplot_full.append(c_condition2)

            pval = poisson_exact(c_condition1, c_condition2)
            if pval < 0.05:
                mark = '*'
            else:
                mark = ''
            marks.append(mark)
            marks.append(mark)

        xjoin = [i for i in range(len(toplot_full))]

        axs[3].bar(xjoin, toplot_full, color=colors_joinbar, width=0.5, linewidth=0)

        start = 0
        bias = ''
        labels = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        for indx, pairs in enumerate(chunks(toplot_full, 2)):
            first, second = pairs[0], pairs[1]
            if (first + second) > start:
                start = first + second
                bias = '{}\t{}\t{}\n'.format(labels[indx], first, second)
        with open('{}/processes/{}/{}.{}.max_components.tsv'.format(outpath, ttype, ttype, s), 'wt') as outfile:
            outfile.write(bias)

        axs[3].set_xticks(xjoin)
        for ix, (mark, val) in enumerate(zip(marks, toplot_full)):
            axs[3].text(xjoin[ix], val+10, mark)

        x = [i for i in range(len(colors))]

        axs[0].bar(x, sig, color=colors, width=0.5, linewidth=0)
        axs[0].set_xticks(x)
        axs[0].set_xticklabels(final_order, rotation=90, fontsize=2)
        axs[0].spines['top'].set_visible(False)

        axs[0].set_xlim(-1, 192)

        colors = []
        colors_mut = ['#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5']
        vector = []
        for lag, lead in zip(sig[0::2], sig[1::2]):
            vector.append(lag-lead)

        for ix, c in enumerate(chunks(vector, 16)):
            colors.extend([colors_mut[ix] for s in c])

        x = [i for i in range(len(colors))]

        axs[1].bar(x, vector, color=colors, width=0.5, linewidth=0)
        axs[1].set_xticks(x)
        axs[1].set_xticklabels(order_plot, rotation=90, fontsize=2)
        axs[1].set_xlim(-1, 96)
        axs[1].spines['top'].set_visible(False)

        vector = []
        colors = []

        sig = sig/np.sum(sig)
        for lag, lead in zip(sig[0::2], sig[1::2]):
            vector.append(lag + lead)

        for ix, c in enumerate(chunks(vector, 16)):
            colors.extend([colors_mut[ix] for s in c])

        x = [i for i in range(len(colors))]
        axs[2].bar(x, vector, color=colors, width=0.75, linewidth=0)
        axs[2].set_xticks(x)
        axs[2].set_xticklabels(order_plot, rotation=90, fontsize=2)
        axs[2].set_xlim(-1, 96)
        axs[2].spines['top'].set_visible(False)

        axs[0].spines['top'].set_visible(False)

        plt.setp([axs[0].get_xticklines(), axs[0].get_yticklines()], color='grey')
        plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')
        plt.setp([axs[2].get_xticklines(), axs[2].get_yticklines()], color='grey')
        plt.setp([axs[3].get_xticklines(), axs[2].get_yticklines()], color='grey')

        axs[0].xaxis.set_ticks_position('none')
        axs[1].xaxis.set_ticks_position('none')
        axs[2].xaxis.set_ticks_position('none')
        axs[3].xaxis.set_ticks_position('none')

        for axis in ['top', 'bottom', 'left', 'right']:
            axs[0].spines[axis].set_linewidth(0.2)
            axs[1].spines[axis].set_linewidth(0.2)
            axs[2].spines[axis].set_linewidth(0.2)
            axs[3].spines[axis].set_linewidth(0.2)
        for indx in [0, 1, 2, 3]:
            axs[indx].xaxis.set_tick_params(pad=0.5)
            axs[indx].yaxis.set_tick_params(pad=0.5, width=0.5)

        plt.tight_layout()

        plt.savefig('{}/processes/{}/{}.{}.png'.format(outpath, ttype, ttype, s), dpi=300, bbox_inches='tight')
        plt.savefig('{}/processes/{}/{}.{}.svg'.format(outpath, ttype, ttype, s))

        plt.close()

        fig, ax = plt.subplots(1, 1, figsize=(1, 1))
        total_max = np.max([np.max(d_clustered[s]), np.max(d_unclustered[s])])

        plt.plot([0, total_max], [0, total_max], lw=1, alpha=0.4)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.scatter(d_clustered[s], d_unclustered[s], s=5, c=colors)
        plt.xlim(0, np.max(d_clustered[s])+1000)
        plt.ylim(0, np.max(d_unclustered[s])+1000)

        plt.xlabel(label1)
        plt.ylabel(label2)

        plt.tight_layout()
        plt.savefig('{}/processes/{}/{}.{}.diagonal.png'.format(outpath, ttype, ttype, s), dpi=300,
                    bbox_inches='tight')
        plt.savefig('{}/processes/{}/{}.{}.diagonal.svg'.format(outpath, ttype, ttype, s))

        plt.close()


# function to plot dbs in bias
def plot_bias_dbs(d_clustered, d_unclustered, outpath, ttype, label1, label2):

    config_params(3)

    final_order = []
    order_plot = order_to_plot_dbs()

    for o in order_plot:
        final_order.append('{}_{}'.format(o, label1))
        final_order.append('{}_{}'.format(o, label2))

    for s in d_clustered.columns:

        sig = d_clustered[s]

        dbs_color = {'AC': '#a6cee3', 'AT': '#1f78b4', 'CC': '#b2df8a', 'CG': '#33a02c', 'CT': '#fb9a99',
                     'GC': '#e3211d', 'TA': '#fdbf6f', 'TC': '#ff7f00', 'TG': '#cab2d6',
                     'TT': '#6a3d9a', }

        fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10, 4))

        vals = sig
        colors = []

        for db in order_plot:
            colors.append(dbs_color[db.split('_')[0]])
            colors.append(dbs_color[db.split('_')[0]])

        x = [i for i in range(len(colors))]

        sig = []
        for lag, lead in zip(d_clustered[s], d_unclustered[s]):
            sig.append(lag)
            sig.append(lead)

        toplot_full = []
        colors_joinbar = []
        marks = []

        start = 0
        for k, color_db in dbs_color.items():
            len_ks = [db for db in order_plot if db.split('_')[0] == k]
            amount = len(len_ks)*2
            c = sig[start:start+amount]
            start = start+amount
            c_condition1 = np.sum(c[::2])
            c_condition2 = np.sum(c[1::2])
            colors_joinbar.append(color_db)
            colors_joinbar.append(color_db)
            toplot_full.append(c_condition1)
            toplot_full.append(c_condition2)

            pval = poisson_exact(c_condition1, c_condition2)
            if pval < 0.05:
                mark = '*'
            else:
                mark = ''
            marks.append(mark)
            marks.append(mark)

        xjoin = [i for i in range(len(toplot_full))]
        axs[3].bar(xjoin, toplot_full, color=colors_joinbar, width=0.5, linewidth=0)
        axs[3].set_xticks(xjoin)
        for ix, (mark, val) in enumerate(zip(marks, toplot_full)):
            axs[3].text(xjoin[ix], val + 10, mark)

        start = 0
        bias = ''
        labels = list(dbs_color.keys())
        for indx, pairs in enumerate(chunks(toplot_full, 2)):
            first, second = pairs[0], pairs[1]
            if (first + second) > start:
                start = first + second
                bias = '{}\t{}\t{}\n'.format(labels[indx], first, second)
        with open('{}/processes/{}/{}.{}.max_components.tsv'.format(outpath, ttype, ttype, s), 'wt') as outfile:
            outfile.write(bias)

        axs[0].bar(x, sig, color=colors, width=0.5, linewidth=0)
        axs[0].set_xticks(x)
        axs[0].set_xticklabels(final_order, rotation=90, fontsize=2)
        axs[0].spines['top'].set_visible(False)

        axs[0].set_xlim(-1, 156)

        colors = []
        vector = []
        for lag, lead in zip(sig[0::2], sig[1::2]):
            vector.append(lag-lead)

        for db in order_plot:
            colors.append(dbs_color[db.split('_')[0]])

        x = [i for i in range(len(colors))]

        axs[1].bar(x, vector, color=colors, width=0.5, linewidth=0)
        axs[1].set_xticks(x)
        axs[1].set_xticklabels(order_plot, rotation=90, fontsize=2)
        axs[1].set_xlim(-1, 78)
        axs[1].spines['top'].set_visible(False)

        vector = []
        colors = []
        sig = sig/np.sum(sig)
        for lag, lead in zip(sig[0::2], sig[1::2]):
            vector.append(lag + lead)

        for db in order_plot:
            colors.append(dbs_color[db.split('_')[0]])

        x = [i for i in range(len(colors))]
        axs[2].bar(x, vector, color=colors, width=0.75, linewidth=0)
        axs[2].set_xticks(x)
        axs[2].set_xticklabels(order_plot, rotation=90, fontsize=2)
        axs[2].set_xlim(-1, 78)

        axs[2].spines['top'].set_visible(False)

        axs[0].spines['top'].set_visible(False)
        axs[0].set_ylabel('NMF counts')
        axs[1].set_ylabel('{} - {}'.format(label1, label2))

        axs[2].set_ylabel('Relative Probability')

        plt.setp([axs[0].get_xticklines(), axs[0].get_yticklines()], color='grey')
        plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')
        plt.setp([axs[2].get_xticklines(), axs[2].get_yticklines()], color='grey')
        plt.setp([axs[3].get_xticklines(), axs[2].get_yticklines()], color='grey')

        axs[0].xaxis.set_ticks_position('none')
        axs[1].xaxis.set_ticks_position('none')
        axs[2].xaxis.set_ticks_position('none')
        axs[3].xaxis.set_ticks_position('none')

        for axis in ['top', 'bottom', 'left', 'right']:
            axs[0].spines[axis].set_linewidth(0.2)
            axs[1].spines[axis].set_linewidth(0.2)
            axs[2].spines[axis].set_linewidth(0.2)
            axs[3].spines[axis].set_linewidth(0.2)
        for indx in [0, 1, 2, 3]:
            axs[indx].xaxis.set_tick_params(pad=0.5)
            axs[indx].yaxis.set_tick_params(pad=0.5, width=0.5)

        plt.tight_layout()

        plt.savefig('{}/processes/{}/{}.{}.png'.format(outpath, ttype, ttype, s), dpi=300, bbox_inches='tight')
        plt.savefig('{}/processes/{}/{}.{}.svg'.format(outpath, ttype, ttype, s))

        plt.close()

        fig, ax = plt.subplots(1, 1, figsize=(1, 1))
        total_max = np.max([np.max(d_clustered[s]), np.max(d_unclustered[s])])

        plt.plot([0, total_max], [0, total_max], lw=1, alpha=0.4)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.scatter(d_clustered[s], d_unclustered[s], s=5, c=colors)
        plt.xlim(0, np.max(d_clustered[s])+10)
        plt.ylim(0, np.max(d_unclustered[s])+10)

        plt.xlabel(label1)
        plt.ylabel(label2)
        plt.tight_layout()

        plt.savefig('{}/processes/{}/{}.{}.diagonal.png'.format(outpath, ttype, ttype, s), dpi=300, bbox_inches='tight')
        plt.savefig('{}/processes/{}/{}.{}.diagonal.svg'.format(outpath, ttype, ttype, s))

        plt.close()


# function to plot indels in bias
def plot_bias_indel(d_clustered, d_unclustered, outpath, ttype, label1, label2):

    config_params(3)

    final_order = []
    order_plot = order_to_plot_indel()

    for o in order_plot:
        final_order.append('{}_{}'.format(o, label1))
        final_order.append('{}_{}'.format(o, label2))

    for s in d_clustered.columns:

        sig = d_clustered[s]

        colors_ind = ['#fdbe6f'] * 6 + ['#ff8002'] * 6 + ['#b0dd8b'] * 6 + ['#36a12e'] * 6 + \
                     ['#fdcab5'] * 6 + ['#fc8a6a'] * 6 + ['#f14432'] * 6 + ['#bc191a'] * 6 + \
                     ['#d0e1f2'] * 6 + ['#94c4df'] * 6 + ['#4a98c9'] * 6 + ['#1764ab'] * 6 + \
                     ['#e1e1ef'] * 1 + ['#b6b6d8'] * 2 + ['#8683bd'] * 3 + ['#62409b'] * 5

        colors_individual = ['#fdbe6f', '#ff8002', '#b0dd8b', '#36a12e', '#fdcab5', '#fc8a6a',
                             '#f14432', '#bc191a', '#d0e1f2', '#94c4df', '#4a98c9', '#1764ab',
                             '#e1e1ef', '#b6b6d8', '#8683bd', '#62409b']

        fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10, 4))

        colors = []

        for col in colors_ind:
            colors.append(col)
            colors.append(col)


        x = [i for i in range(len(colors))]

        sig = []
        for lag, lead in zip(d_clustered[s], d_unclustered[s]):
            sig.append(lag)
            sig.append(lead)

        axs[0].bar(x, sig, color=colors, width=0.8, linewidth=0)
        axs[0].set_xticks(x)
        axs[0].set_xticklabels(final_order, rotation=90, fontsize=2)
        axs[0].spines['top'].set_visible(False)

        axs[0].set_xlim(-1, 172)

        toplot_full = []
        colors_joinbar = []
        marks = []

        start = 0

        sizes = [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 2, 3, 5]

        for ix, len_size in enumerate(sizes):

            c = sig[start:start+len_size * 2]
            start = start+len_size * 2
            c_condition1 = np.sum(c[::2])
            c_condition2 = np.sum(c[1::2])
            colors_joinbar.append(colors_individual[ix])
            colors_joinbar.append(colors_individual[ix])
            toplot_full.append(c_condition1)
            toplot_full.append(c_condition2)

            pval = poisson_exact(c_condition1, c_condition2)
            if pval < 0.05:
                mark = '*'
            else:
                mark = ''
            marks.append(mark)
            marks.append(mark)

        xjoin = [i for i in range(len(toplot_full))]
        axs[3].bar(xjoin, toplot_full, color=colors_joinbar, width=0.5, linewidth=0)
        axs[3].set_xticks(xjoin)
        for ix, (mark, val) in enumerate(zip(marks, toplot_full)):
            axs[3].text(xjoin[ix], val + 10, mark)

        vector = []
        for lag, lead in zip(sig[0::2], sig[1::2]):
            vector.append(lag-lead)

        x = [i for i in range(len(colors_ind))]

        axs[1].bar(x, vector, color=colors_ind, width=0.5, linewidth=0)
        axs[1].set_xticks(x)
        axs[1].set_xticklabels(order_plot, rotation=90, fontsize=2)
        axs[1].set_xlim(-1, 86)
        axs[1].spines['top'].set_visible(False)

        vector = []

        sig = sig/np.sum(sig)
        for lag, lead in zip(sig[0::2], sig[1::2]):
            vector.append(lag + lead)

        x = [i for i in range(len(colors_ind))]
        axs[2].bar(x, vector, color=colors_ind, width=0.75, linewidth=0)
        axs[2].set_xticks(x)
        axs[2].set_xticklabels(order_plot, rotation=90, fontsize=2)
        axs[2].set_xlim(-1, 86)

        axs[2].spines['top'].set_visible(False)

        axs[0].spines['top'].set_visible(False)
        axs[0].set_ylabel('NMF counts')
        axs[1].set_ylabel('{} - {}'.format(label1, label2))

        axs[2].set_ylabel('Relative Probability')

        plt.setp([axs[0].get_xticklines(), axs[0].get_yticklines()], color='grey')
        plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')
        plt.setp([axs[2].get_xticklines(), axs[2].get_yticklines()], color='grey')
        plt.setp([axs[3].get_xticklines(), axs[2].get_yticklines()], color='grey')

        axs[0].xaxis.set_ticks_position('none')
        axs[1].xaxis.set_ticks_position('none')
        axs[2].xaxis.set_ticks_position('none')
        axs[3].xaxis.set_ticks_position('none')

        for axis in ['top', 'bottom', 'left', 'right']:
            axs[0].spines[axis].set_linewidth(0.2)
            axs[1].spines[axis].set_linewidth(0.2)
            axs[2].spines[axis].set_linewidth(0.2)
            axs[3].spines[axis].set_linewidth(0.2)
        for indx in [0, 1, 2, 3]:
            axs[indx].xaxis.set_tick_params(pad=0.5)
            axs[indx].yaxis.set_tick_params(pad=0.5, width=0.5)

        plt.tight_layout()

        plt.savefig('{}/processes/{}/{}.{}.png'.format(outpath, ttype, ttype, s), dpi=300, bbox_inches='tight')
        plt.savefig('{}/processes/{}/{}.{}.svg'.format(outpath, ttype, ttype, s))

        plt.close()

        fig, ax = plt.subplots(1, 1, figsize=(1, 1))
        total_max = np.max([np.max(d_clustered[s]), np.max(d_unclustered[s])])

        plt.plot([0, total_max], [0, total_max], lw=1, alpha=0.4)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.scatter(d_clustered[s], d_unclustered[s], s=5, c=colors_ind)
        plt.xlim(0, np.max(d_clustered[s])+10)
        plt.ylim(0, np.max(d_unclustered[s])+10)

        plt.xlabel(label1)
        plt.ylabel(label2)

        plt.tight_layout()
        plt.savefig('{}/processes/{}/{}.{}.diagonal.png'.format(outpath, ttype, ttype, s), dpi=300,
                    bbox_inches='tight')
        plt.savefig('{}/processes/{}/{}.{}.diagonal.svg'.format(outpath, ttype, ttype, s))

        plt.close()
