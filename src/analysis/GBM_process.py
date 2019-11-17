# Import modules
import os
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from utils.plots_utils import config_params


def get_tmz_deconstructsigs():
    decons_file = "data/GBM_Wang/gbm_muts.to_deconstruct.out.tsv"
    dfdec = pd.read_csv(decons_file, sep='\t', index_col=0)

    total_rows = []
    for sample, row in dfdec.iterrows():
        out = row['mutation_count']*row
        total_rows.append(out)
    extracted_sig = pd.DataFrame(total_rows)
    return extracted_sig.T


def read_clinical_data():

    clinic = 'data/GBM_Wang/ng.3590-S12.xlsx'
    dfclin = pd.read_excel(clinic)

    good_treatment = set()
    for i, row in dfclin.iterrows():
        if ('Stupp' in row['Non-Surgical Treatment']) or ('TMZ' in row['Non-Surgical Treatment']):
            good_treatment.add(row['Non-Surgical Treatment'])

    # MGMT in recurrence
    total_MGMT = set(dfclin[dfclin['MGMT_Recurrence'] == 1]['#ID'].tolist())

    wanted = ['CCRT', 'None', 'RT, Procarbazine/Vincristine/Lomustine']
    nottreated = dfclin[dfclin['Non-Surgical Treatment'].isin(wanted)]['#ID']
    total_treated = dfclin[
        (dfclin['Non-Surgical Treatment'].isin(good_treatment)) &
        (dfclin['Type'] == 'Secondary GBM')
        ]['#ID']
    return total_MGMT, nottreated, total_treated


def get_MMR_def():

    rep = 'data/repair_genes/repair_publication.tsv'
    clinic = 'data/GBM_Wang/ng.3590-S4.xlsx'
    dfmuts = pd.read_excel(clinic, sep='\t')
    df_rep = pd.read_csv(rep, sep='\t', names=['gene', 'pathway'])
    mmr_genes = df_rep["gene"].tolist()
    samples_mmr = set(
        dfmuts[
            (dfmuts['Gene_Name'].isin(mmr_genes)) &
            (dfmuts['Effect'] != 'synonymous_variant')
            ]['CaseID'].tolist()
    )
    return samples_mmr


def GBM_tzm_plot(total_MGMT, nottreated, total_treated, samples_mmr, fitted):

    sig_temo = 'SBS11'
    config_params(6.5)

    # first we select MMR defitient
    treated_MMR = [i for i in total_treated if i in samples_mmr]

    # then we select MGMT but nor MMR deff
    treated_MGMT = [i for i in total_treated if ((i not in treated_MMR) & (i in total_MGMT))]

    # then treated without the others
    treated_no_alt = [i for i in total_treated if (i not in treated_MMR) & (i not in treated_MGMT)]

    fig, ax = plt.subplots(1, 1, figsize=(1.25, 1.5))

    sns.boxplot(
        data=[
            fitted[[i for i in nottreated if i in fitted.columns]].loc[sig_temo],
            fitted[[i for i in treated_no_alt if i in fitted.columns]].loc[sig_temo],
            fitted[[i for i in treated_MGMT if i in fitted.columns]].loc[sig_temo],
            fitted[[i for i in treated_MMR if i in fitted.columns]].loc[sig_temo]
        ],
        color='#ecececff',
        linewidth=1,
        showfliers=False,
        ax=ax)

    ax.scatter(
        [0 + np.random.uniform(-0.2, 0.2, 1)[0] for i in fitted[
            [i for i in nottreated if i in fitted.columns]
        ].loc[sig_temo]],
        fitted[[i for i in nottreated if i in fitted.columns]].loc[sig_temo],
        s=15,
        color='#800080ff',
        edgecolor='black',
        lw=0.5
    )

    ax.scatter(
        [1 + np.random.uniform(-0.2, 0.2, 1)[0] for i in fitted[
            [i for i in treated_no_alt if i in fitted.columns]
        ].loc[sig_temo]],
        fitted[[i for i in treated_no_alt if i in fitted.columns]].loc[sig_temo],
        s=15,
        color='#800080ff',
        edgecolor='black',
        lw=0.5
    )

    ax.scatter(
        [2 + np.random.uniform(-0.2, 0.2, 1)[0] for i in fitted[
            [i for i in treated_MGMT if i in fitted.columns]
        ].loc[sig_temo]],
        fitted[[i for i in treated_MGMT if i in fitted.columns]].loc[sig_temo],
        s=15,
        color='#800080ff',
        edgecolor='black',
        lw=0.5
    )

    ax.scatter(
        [3 + np.random.uniform(-0.1, 0.1, 1)[0] for i in fitted[
            [i for i in treated_MMR if i in fitted.columns]
        ].loc[sig_temo]],
        fitted[[i for i in treated_MMR if i in fitted.columns]].loc[sig_temo],
        s=15,
        color='#800080ff',
        edgecolor='black',
        lw=0.5
    )

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.xticks(
        [0, 1, 2, 3],
        [
            'TMZ untreated (n ={})'.format(len(nottreated)),
            'MGMT-notmet MMR-prof (n = {})'.format(len(treated_no_alt)),
            'MGMT-met MMR-prof (n = {})'.format(len(treated_MGMT)),
            'MMR-def(n = {})'.format(len(treated_MMR)),
        ],
        rotation=90
    )

    plt.title('GBM-cohort\n(Wang et al, 2016)')
    plt.ylabel('TMZ related SBS')
    plt.savefig('figures/GBM_tmz.svg')
    plt.close()


def run():
    total_MGMT, nottreated, total_treated = read_clinical_data()
    fitted = get_tmz_deconstructsigs()
    samples_mmr = get_MMR_def()
    GBM_tzm_plot(total_MGMT, nottreated, total_treated, samples_mmr, fitted)


if __name__ == '__main__':
    run()
