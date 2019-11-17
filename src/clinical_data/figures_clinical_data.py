# Import modules
import sys
import os
from collections import defaultdict
import gzip
import pickle

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from tqdm import tqdm
import holoviews as hv

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)
scripts_path = os.path.abspath(os.path.join(__file__, '..', '..', 'clinical_data/'))
sys.path.insert(0, scripts_path)
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf
from utils.metadata import return_metadata
from utils.colors_ttype import return_colors
from utils.plots_utils import config_params


# Set the random seed for reproducibility
np.random.seed(12345)


def create_matrix_treatments_plot():

    treated = pickle.load(gzip.open(Conf['treatment_FDA_drug']))
    dic_primary, _ = return_metadata()
    matrix_treated = defaultdict(lambda: defaultdict(int))

    forbidden = ['CUP', 'Eye', 'Double-primary', 'Unknown']

    for sample, t in dic_primary.items():
        if t not in forbidden:
            for k, d in treated.items():
                if sample in d[t]['YES']:
                    matrix_treated[sample][k] = 1
                elif sample in d[t]['NO']:
                    matrix_treated[sample][k] = 0

    out = pd.DataFrame.from_dict(matrix_treated, orient='index')
    out = out.dropna()
    dic_t = defaultdict(list)
    for i, row in out.iterrows():
        dic_t[dic_primary[i]].append(i)

    return out, dic_t


def sankey_plot_main():

    config_params(font_size=4)

    hv.extension('matplotlib')
    hv.output(fig='svg')

    forbidden = ['RADIATION', 'Miscellanious', 'Unknown', 'TopoII', 'TOPOII']
    out, dic_t = create_matrix_treatments_plot()
    order_ttypes = [
        'Breast',
        'Colon-Rectum',
        'Prostate',
        'Lung',
        'Skin',
        'Bone-Soft-tissue',
        'Ovary',
        'Esophagus',
        'Urinary-tract',
        'NET',
        'Kidney',
        'Nervous-system',
        'Biliary',
        'Pancreas',
        'Unknown',
        'Uterus',
        'Head-and-neck',
        'Liver',
        'Stomach',
        'Mesothelioma',
    ]

    all_rows = []
    for ttype in order_ttypes:
        samples = dic_t[ttype]
        subs = out.loc[samples]
        for col in subs:
            if col not in forbidden:
                all_rows.append((ttype, col, int(subs[col].sum())))

    matrix_df = pd.DataFrame(all_rows)
    matrix_df.columns = ['target', 'source',  'value']
    matrix_df = matrix_df[(matrix_df['target'] != 'Unknown')]
    matrix_df = matrix_df.fillna(0)
    matrix_df['value'] = matrix_df['value'].astype(int)

    good_source = set()
    for source, data in matrix_df.groupby(by='source'):
        tot = data['value'].sum()
        if tot > 30:
            if source != 'Unknown':
                good_source.add(source)
    matrix_df = matrix_df[matrix_df['source'].isin(good_source)]
    out = hv.Sankey(matrix_df.sort_values(by='source', ascending=True, ), label='').opts(
        label_position='left', edge_color='target', node_color='index', cmap='Set1')  # color=total_colors)

    fig = hv.render(out)
    fig.set_figwidth(10)
    fig.savefig('figures/2A.svg')
    fig.savefig('figures/2A.png', dpi=600)


def plot_distribution_tissue_origin():

    # PLOT DISTRIBUTION NUMBER COHORTS
    file_metadata = Conf['path_metadata']

    colors_ttype = return_colors()

    # read Hartwig metadata
    pat = pd.read_csv(file_metadata, sep='\t')

    pat['primaryTumorLocation'] = pat['primaryTumorLocation'].replace('Bone/soft tissue', 'Bone/Soft tissue')

    # fix primary location
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation'].apply(
        lambda x: str(x).replace(' ', '-').replace('/', '-')
    )
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('Head-and-Neck', 'Head-and-neck')
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('nan', 'Unknown')
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('CUP', 'Unknown')
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('Net', 'NET')

    dic_tumor = pat['primaryTumorLocation_fixed'].value_counts().to_dict()

    sorted_keys = sorted(dic_tumor, key=dic_tumor.get, reverse=True)

    config_params(font_size=5.5)
    fig, ax = plt.subplots(1, 1, figsize=(2, 2.35))
    labels = []
    count = 0
    for ix, k in enumerate(sorted_keys[::-1]):
        if dic_tumor[k] > 20:
            if k in colors_ttype:
                count += 1
                ax.barh(count, dic_tumor[k], color=colors_ttype[k])
                ax.text(dic_tumor[k] + 40, count - 0.1, dic_tumor[k], verticalalignment='center')
                labels.append(k)

    plt.yticks([i + 1 for i in range(0, len(labels))], labels, rotation=0)

    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.xlabel('Number of samples')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    plt.tight_layout()
    ax.xaxis.set_tick_params(width=0.45)
    ax.yaxis.set_tick_params(width=0.45)
    plt.savefig('figures/1B.svg')
    plt.savefig('figures/1B.png', dpi=600)

    plt.show()


def specific_treatments_bar(drug):

    config_params(5)
    dplot_treat = defaultdict(int)
    typeDrug = pickle.load(gzip.open('data/clinical_data/hartwig_typeDrug.pckl.gz'))
    samples_tract_specific_FDA = pickle.load(gzip.open(Conf['treatment_specific_drug']))

    for k in typeDrug[drug].keys():
        number_treated = samples_tract_specific_FDA[k]['Pan']['YES']
        dplot_treat[k] = len(number_treated)

    sorted_keys = sorted(dplot_treat, key=dplot_treat.get, reverse=True)

    fig, ax = plt.subplots(1, 1, figsize=(0.3, 1))
    bottom = 0
    for treat in sorted_keys:
        ax.bar(0, dplot_treat[treat], width=1, bottom=bottom)
        ax.text(-0.1, dplot_treat[treat] + bottom - 150, str(dplot_treat[treat]))
        bottom += dplot_treat[treat]

    plt.xticks([0], fontsize=5)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig('figures/2A_{}_barplot.svg'.format(drug.lower().replace(' ', '_')))
    plt.savefig('figures/2A_{}_barplot.png'.format(drug.lower().replace(' ', '_')), dpi=600)

    plt.close()


def classic_mutual_exclusivity_visualization(regression_table, sorted_treatments):

    df = regression_table.copy()
    df = df[sorted_treatments]
    X = df.values
    # for each gene give weights proportional to the number of mutations there
    weights = [1 / (i + 1) for i, c in enumerate(df.columns)]
    Y = pdist(X, metric='hamming', w=np.array(weights))
    linkage = hierarchy.linkage(Y, method='ward')
    g = sns.clustermap(X.T, col_linkage=linkage, col_cluster=True, row_cluster=False, figsize=(20, 10), )
    g.cax.set_visible(False)
    plt.close()

    return g


def plot_heatmap_treatment(file):

    outf = os.path.basename(file).split('.')[0]
    dic_primary_full, _ = return_metadata()
    color_ttype = return_colors()
    total_s = defaultdict(list)
    total_count = defaultdict(int)

    for sample, t in dic_primary_full.items():
        total_s[t].append(sample)
        total_count[t] += 1
    sorted_ttyps = sorted(total_count, key=total_count.get, reverse=True)

    treated = pickle.load(gzip.open(file))
    forbidden = ['RADIATION', 'TOPOII']
    matrix_treated = defaultdict(lambda: defaultdict(int))

    for sample, t in dic_primary_full.items():
        for k, d in treated.items():
            if k not in forbidden:
                if sample in d[t]['YES']:
                    matrix_treated[sample][k] = 1
                elif sample in d[t]['NO']:
                    matrix_treated[sample][k] = 0

    d_treatments = defaultdict(int)
    for k, d in treated.items():
        if k not in forbidden:
            for ttype, l in d.items():
                d_treatments[k] += len(l['YES'])

    sorted_treatments = sorted(d_treatments, key=d_treatments.get, reverse=True)
    out = pd.DataFrame.from_dict(matrix_treated, orient='index')
    out['TTYPE'] = [dic_primary_full[t] for t in out.index.tolist()]
    out['sum'] = out.sum(axis=1)
    out = out[out['sum'] > 0].drop('sum', axis=1)

    forbidden = ['Double-primary']
    order_sample_plot = []
    order = []
    dic_len = defaultdict(int)
    for ttype in tqdm(sorted_ttyps):
        if ttype not in forbidden:
            subs = out[out['TTYPE'] == ttype]
            mat = subs.drop('TTYPE', axis=1).dropna()[sorted_treatments[:30]]
            if len(mat) > 1:
                n = classic_mutual_exclusivity_visualization(mat, sorted_treatments[:30])
                new_order = n.dendrogram_col.reordered_ind
                sample_list = mat.reset_index().loc[new_order]['index'].tolist()
                order_sample_plot.extend(sample_list)
                order.append(ttype)
                dic_len[ttype] = len(sample_list)

    new_cmap = LinearSegmentedColormap.from_list("", ["lightgrey", "grey", "darkred"])
    concat = out.loc[order_sample_plot].drop('TTYPE', axis=1)
    concat = concat[sorted_treatments[:20]]

    if 'specific' in outf:
        new_cols = [s.lower() for s in concat.columns]
        concat.columns = new_cols

    config_params(2)
    fig, ax = plt.subplots(1, 2, figsize=(1, 3), gridspec_kw={'width_ratios': [1, 27]})
    ax2 = sns.heatmap(
        concat, cmap=new_cmap, yticklabels=False, ax=ax[1], cbar=False
    )
    ax[1].xaxis.set_ticks_position('top')
    bot = 0
    for t in order[::-1]:
        ax[0].bar(0, dic_len[t], bottom=bot, color=color_ttype[t])
        bot += dic_len[t]
    ax[0].set_ylim(0, bot)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)
    ax[0].spines['left'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].get_yaxis().set_visible(False)
    ax[0].get_xaxis().set_visible(False)

    plt.xticks(rotation=90)
    plt.savefig('figures/EDF1_{}.png'.format(outf), dpi=600)
    plt.close()


if __name__ == '__main__':

    # Figure 1B
    plot_distribution_tissue_origin()

    # Figure 2A
    sankey_plot_main()
    specific_treatments_bar('Platinum-based Drug')
    specific_treatments_bar('Nucleoside Metabolic Inhibitor')

    # Supp Figures
    plot_heatmap_treatment(Conf['treatment_FDA_drug'])
    plot_heatmap_treatment(Conf['treatment_specific_drug'])
