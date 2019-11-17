# Import modules
import os
import sys
from glob import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from utils.plots_utils import config_params


def get_results(method_extraction, d_sig):

    all_regs = os.path.join(
        'data', 'regression_treatments', method_extraction,
        'merged', 'regression_results_treatments_specific_*.tsv'
    )

    size_type = {
        'snvs': "o",
        'dbs': 'v',
        'indels': "P",
    }

    color_treatment = {
        'oxaliplatin': '#ffb380ff',
        'carboplatin': '#ff6600ff',
        'cisplatin': '#ffe6d5ff',
        'capecitabine': '#2c89a0ff',
        'temozolomide': 'black'
    }

    toconcat = []
    for f in glob(all_regs):
        df = pd.read_csv(f, sep='\t')
        type_var = os.path.basename(f).split('_')[-1].split('.')[0]
        type_var_final = type_var
        if type_var == 'snv':
            type_var_final = 'snvs'
        df['TYPE'] = type_var_final
        df['sig_unique'] = df['signature'] + "_" + df['TYPE']
        df['format_type'] = df['TYPE'].map(size_type)
        df['pvals'] = df['pvals'].replace(0, 0.0001)
        df['logp'] = df['pvals'].apply(lambda x: -np.log10(x))
        toconcat.append(df)

    all_regressions = pd.concat(toconcat)
    significant = all_regressions[(all_regressions['effect_size'] > 2) & (all_regressions['pvals'] < 0.001)]

    nonsignificant = all_regressions[(all_regressions['effect_size'] <= 2) | (all_regressions['pvals'] >= 0.001)]

    all_class = []
    for i, row in significant.iterrows():
        class_rel = 'VETTED'
        if row['treatment'] in d_sig:
            if row['sig_unique'] in d_sig[row['treatment']]:
                class_rel = 'OK'
        all_class.append(class_rel)

    significant['VETTING'] = all_class
    significant['color'] = significant.apply(
        lambda x: color_treatment[x['treatment']] if x['VETTING'] == 'OK' else '#99d8c9', axis=1
    )
    significant_not_vetted = significant[significant['VETTING'] == 'OK']
    significant_vetted = significant[significant['VETTING'] != 'OK']

    nonsignificant['VETTING'] = 'NON_SIG'
    nonsignificant['color'] = 'grey'
    do_volcano(significant_vetted, nonsignificant, significant_not_vetted, method_extraction)


def do_volcano(significant_vetted, nonsignificant, significant_not_vetted, method_extraction):

    config_params()
    fig, ax = plt.subplots(1, 1, figsize=(4.2, 4))
    size = 3
    plt.hlines(-np.log10(0.001), 0, 10, linestyles="dashed", color='red', alpha=0.3)
    plt.vlines(2, 0, 7, linestyles="dashed", color='red', alpha=0.3)

    for i, row in significant_vetted.iterrows():
        plt.scatter(row['effect_size'], row['logp'],
                    c='#99d8c9', s=size, marker=row['format_type'])

    for i, row in nonsignificant.iterrows():
        plt.scatter(row['effect_size'], row['logp'], c='grey',
                    s=size, marker=row['format_type'])

    for i, row in significant_not_vetted.iterrows():
        plt.scatter(row['effect_size'], row['logp'],
                    c=row['color'], s=25, edgecolor='black', linewidths=0.7,
                    marker=row['format_type'])

    plt.xlim(0, 16)
    plt.ylim(0, 3.1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.ylabel('statistical significance\n(-log10 pvalue)')
    plt.xlabel('treated-untreated fold change')
    plt.tight_layout()
    plt.savefig('figures/{}/volcano.svg'.format(method_extraction))
    plt.savefig('figures/{}/volcano.png'.format(method_extraction), dpi=600)

    plt.show()


def do_plot():

    d_sig_signature_analyzer = {
        'cisplatin': ['21_SBS31_0.953955_1_snvs', '14_1_snvs', '9_1_dbs', '3_1_dbs'],
        'carboplatin': ['21_SBS31_0.953955_1_snvs', '25_1_snvs', '9_1_dbs', '3_1_dbs'],
        'capecitabine': ['31_SBS17b_0.968799_1_snvs'],
        'oxaliplatin': ['14_1_snvs', '37_1_snvs', '9_1_dbs', '3_1_dbs'],
    }
    d_sig_sigprofiler = {
        'cisplatin': ['1_SBS31_0.968153_0.98_snvs', '5_DBS5_0.944431_1.0_dbs'],
        'carboplatin': ['5_DBS5_0.944431_1.0_dbs'],
        'oxaliplatin': ['20_0.92_snvs', '5_DBS5_0.944431_1.0_dbs'],
        'capecitabine': ['19_SBS17b_0.961548_0.99_snvs'],
    }

    get_results("SigProfiler", d_sig_sigprofiler)
    get_results("SignatureAnalyzer", d_sig_signature_analyzer)


if __name__ == '__main__':
    do_plot()
