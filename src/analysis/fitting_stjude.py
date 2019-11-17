# Import modules
import os
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

sys.exit()

from utils.plots_utils import config_params


def do_plot(list_val, norm=False):

    fig, axs = plt.subplots(1, 1, figsize=(0.1, 1.5))
    toplot = 0
    config_params()

    if norm is False:
        axs.set_xlabel('{} ({})'.format('Pediatric', len(list_val)), rotation=90)
        axs.set_ylim(10, 200000)
        axs.set_yscale('log')

    range_nums = [num / len(list_val) for num in range(len(list_val))]

    axs.hlines(
        np.median(list_val),
        np.median(range_nums) - 0.25, np.median(range_nums) + 0.45,
        color='grey', alpha=0.4
    )
    axs.scatter(
        range_nums,
        sorted(list_val), s=2,
        linewidth=0.2, color='#006400ff'
    )

    axs.set_xlim(-0.2, 1.1)
    axs.spines['bottom'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.set_xticklabels([])

    if toplot > 0:
        axs.spines['left'].set_visible(False)
    else:
        axs.set_ylabel('Number of SBS\nassociated to treatments')

    lab = 'full'
    if norm is True:
        axs.set_ylim(0, 0.8)
        lab = 'norm'

    axs.xaxis.set_label_position('top')
    axs.xaxis.set_ticks_position('none')
    axs.yaxis.set_ticks_position('none')

    plt.savefig('figures/stjude_{}.svg'.format(lab))


def process_stjude_fitting():

    df = pd.read_csv('data/STJUDE/format/STJUDE.to_deconstruct.out.tsv', sep='\t')
    df.set_index('sample_id', inplace=True)

    list_val = []
    list_val_norm = []

    for i, row in df.iterrows():
        sig_cis = (row['SBS31'] + row['SBS35']) * row['mutation_count']
        if sig_cis > 0:
            list_val.append(sig_cis)
            list_val_norm.append((row['SBS31'] + row['SBS35']))

    # do plots
    do_plot(list_val, norm=False)
    do_plot(list_val_norm, norm=True)


if __name__ == '__main__':
    process_stjude_fitting()
