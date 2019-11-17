import pandas as pd
import scipy.spatial.distance as spdist
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib as mpl

# config for matplotlib
def config_params(font_size=7):
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'

# split into even chunks
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

def order_muts():

    order = []
    first = ['A', 'C', 'G', 'T']
    pyr = ['C', 'T']
    for p in pyr:
        for mut in first:
            if mut!=p:
                for f in first:
                    for f2 in first:
                        comb = '{}[{}>{}]{}'.format(f, p, mut, f2)
                        order.append(comb)
    return order


def plot_signatures(sig):

    config_params(4.5)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(3.2, 1),
                            gridspec_kw={'height_ratios': [1, 9]})
    order_plot = order_muts()

    vals = []
    colors = []
    colors_mut = ['#1ebff0', '#050708', '#e62725',
                  '#cbcacb', '#a1cf64', '#edc8c5']
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

    axs[1].axhline(y=0.05, xmin=-1, xmax=96, lw = 0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.1, xmin=-1, xmax=96, lw = 0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.15, xmin=-1, xmax=96, lw = 0.6, color='grey', alpha=0.2)

    axs[1].bar(x, vals, color=colors, width=0.8, linewidth=0, align='center')
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(['{}{}{}'.format(a[0], a[2], a[-1]) for a in order_plot],
                           verticalalignment="center", ha='center', rotation=90, fontsize=3,
                           color='grey')

    plt.tight_layout()
    plt.xlim(-1, 96)

    axs[1].spines['top'].set_visible(False)
    axs[1].set_ylabel('Relative Probability', fontsize = 4.5)
    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')
    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=0.5)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)


    plt.tick_params(axis='both', which='both', bottom=False, left = False)
    plt.show()


def get_cosine_similarity():

    context = ["C.A.in.ACA","C.A.in.ACC","C.A.in.ACG","C.A.in.ACT","C.A.in.CCA",
               "C.A.in.CCC","C.A.in.CCG","C.A.in.CCT","C.A.in.GCA","C.A.in.GCC",
               "C.A.in.GCG","C.A.in.GCT","C.A.in.TCA","C.A.in.TCC","C.A.in.TCG",
               "C.A.in.TCT","C.G.in.ACA","C.G.in.ACC","C.G.in.ACG","C.G.in.ACT",
               "C.G.in.CCA","C.G.in.CCC","C.G.in.CCG","C.G.in.CCT","C.G.in.GCA",
               "C.G.in.GCC","C.G.in.GCG","C.G.in.GCT","C.G.in.TCA","C.G.in.TCC",
               "C.G.in.TCG","C.G.in.TCT","C.T.in.ACA","C.T.in.ACC","C.T.in.ACG",
               "C.T.in.ACT","C.T.in.CCA","C.T.in.CCC","C.T.in.CCG","C.T.in.CCT",
               "C.T.in.GCA","C.T.in.GCC","C.T.in.GCG","C.T.in.GCT","C.T.in.TCA",
               "C.T.in.TCC","C.T.in.TCG","C.T.in.TCT","T.A.in.ATA","T.A.in.ATC",
               "T.A.in.ATG","T.A.in.ATT","T.A.in.CTA","T.A.in.CTC","T.A.in.CTG",
               "T.A.in.CTT","T.A.in.GTA","T.A.in.GTC","T.A.in.GTG","T.A.in.GTT",
               "T.A.in.TTA","T.A.in.TTC","T.A.in.TTG","T.A.in.TTT","T.C.in.ATA",
               "T.C.in.ATC","T.C.in.ATG","T.C.in.ATT","T.C.in.CTA","T.C.in.CTC",
               "T.C.in.CTG","T.C.in.CTT","T.C.in.GTA","T.C.in.GTC","T.C.in.GTG",
               "T.C.in.GTT","T.C.in.TTA","T.C.in.TTC","T.C.in.TTG","T.C.in.TTT",
               "T.G.in.ATA","T.G.in.ATC","T.G.in.ATG","T.G.in.ATT","T.G.in.CTA",
               "T.G.in.CTC","T.G.in.CTG","T.G.in.CTT","T.G.in.GTA","T.G.in.GTC",
               "T.G.in.GTG","T.G.in.GTT","T.G.in.TTA","T.G.in.TTC","T.G.in.TTG","T.G.in.TTT"]

    df_processes_HDP = pd.read_csv('data/hartwig/signatures/HDP/processes_method_HDP.tsv', sep ='\t', decimal=",")
    df_processes_signatureanalyzer = pd.read_csv('data/hartwig/signatures/extraction/results/SignatureAnalyzer/snvs/processes/Pan_full/Pan_full.processes.tsv', sep ='\t',
                                                index_col = 0)

    context_modified = ['{}[{}>{}]{}'.format(x[-3], x[-2], x[2], x[-1]) for x in context]
    index_modified = ['HDP_{}'.format(indx) for indx in df_processes_HDP.index.tolist()]
    df_processes_HDP.columns =  context_modified
    df_processes_HDP.index = index_modified
    df_processes_HDP = df_processes_HDP.T

    df_processes_signatureanalyzer = df_processes_signatureanalyzer.loc[context_modified]
    d_sigan = defaultdict(list)
    for col in df_processes_signatureanalyzer.columns:
        d_sigan[col] = df_processes_signatureanalyzer[col].tolist()


    # find similar signatures to those reported in SigProfiler
    cos_sim = defaultdict(dict)
    for ix, col in enumerate(df_processes_HDP.columns):
        vec1 = df_processes_HDP[col].tolist()
        for s, vec2 in d_sigan.items():
            c = 1 - round(spdist.cosine(vec1, vec2), 7)
            cos_sim[col][s] = c

    # select those with higher similarity per each of the signatures
    cos_df = pd.DataFrame(cos_sim)
    index_max = cos_df.idxmax()
    vals = cos_df.max()
    similar = pd.DataFrame(list(zip(index_max, vals)))

    similar.index = cos_df.columns

    return df_processes_HDP, similar


def multiply_col(df, dict_number_muts):
    col_name = df.name
    return df*dict_number_muts[col_name]

# get full exposures to run the regression afterwards
def get_full_exposures():

    df_exposures_HDP = pd.read_csv('data/hartwig/signatures/HDP/exp_method_HDP.tsv', sep ='\t', decimal=",")

    df_original = pd.read_csv('data/hartwig/signatures/matrix/Colon-Rectum.HDP.txt', sep ='\t', decimal=",",
                             index_col = 0)

    dict_number_muts = df_original.sum(axis = 1).to_dict()
    list_samples = df_original.index.tolist()

    df_exposures_HDP.columns = list_samples

    df_exposures_HDP_full = df_exposures_HDP.apply(multiply_col, args = (dict_number_muts, ))
    df_exposures_HDP_full.reset_index(drop = True, inplace = True)
    df_exposures_HDP_full.index = ['HDP_{}'.format(i) for i in range(1, len(df_exposures_HDP_full)+1)]

    df_exposures_HDP_full.to_csv('data/hartwig/signatures/HDP/exp_method_HDP_full.tsv', sep ='\t', index = True, header = True)


# evaluate cosinus similarity
df_processes_HDP, similar = get_cosine_similarity()

print(similar)

sig = df_processes_HDP['HDP_10'].tolist()
plot_signatures(sig)

sig = df_processes_HDP['HDP_8'].tolist()
plot_signatures(sig)

get_full_exposures()