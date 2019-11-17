# Import modules
import os
import sys
import gzip
from glob import glob
from collections import defaultdict

import pandas as pd
import pathlib
import scipy.spatial.distance as spdist
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from tqdm import tqdm
import dill as pickle

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
scripts_path = os.path.abspath(os.path.join(__file__, '..', '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf
from utils.orders import order_to_plot_indel, order_to_plot_dbs, order_to_plot_snvs, \
    snvs_signature_analyzer, order_sigprofiler_snvs
from utils.signature_plots import plot_snvs, plot_dbs, plot_indel, config_params, \
    plot_bias_indel, plot_bias_dbs, plot_bias_snv


# read  canonical SBS signatures from PCAWG
def pcawg_canonical_snvs():
    pcawg_sbs_file = Conf['SBS_PCAWG']
    pcawg_snvs = pd.read_csv(pcawg_sbs_file)
    pcawg_snvs.index = pcawg_snvs.apply(lambda x: '{}[{}>{}]{}'.format(
        x['SubType'][0],
        x['SubType'][1],
        x['Type'][-1],
        x['SubType'][2]), axis=1
                                        )

    pcawg_snvs.sort_index(inplace=True)
    pcawg_snvs.drop(['SubType', 'Type'], axis=1, inplace=True)
    # load them in order
    d_pcawg = defaultdict(list)

    for col in pcawg_snvs.columns:
        d_pcawg[col] = pcawg_snvs[col].tolist()

    return d_pcawg, pcawg_snvs.index.tolist()


# read canonical DBS signatures from PCAWG
def pcawg_canonical_dbs():
    pcawg_dbs_file = Conf['DBS_PCAWG']
    pcawg_dbs = pd.read_csv(pcawg_dbs_file)
    pcawg_dbs.sort_values(by='Mutation Type', inplace=True)
    pcawg_dbs.index = pcawg_dbs['Mutation Type']
    pcawg_dbs.drop(['Mutation Type'], axis=1, inplace=True)
    d_pcawg = defaultdict(list)

    for col in pcawg_dbs.columns:
        d_pcawg[col] = pcawg_dbs[col].tolist()

    return d_pcawg, pcawg_dbs.index.tolist()


# read canonical indel signatures from PCAWG
def pcawg_canonical_indels():
    pcawg_indels_file = Conf['ID_PCAWG']
    pcawg_indels = pd.read_csv(pcawg_indels_file)
    pcawg_indels.sort_values(by='Mutation Type', inplace=True)
    pcawg_indels.index = pcawg_indels['Mutation Type']
    pcawg_indels.drop(['Mutation Type'], axis=1, inplace=True)
    d_pcawg = defaultdict(list)

    for col in pcawg_indels.columns:
        d_pcawg[col] = pcawg_indels[col].tolist()

    return d_pcawg, pcawg_indels.index.tolist()


# compact channels
def return_1536_to_96_channels(pro_df):
    Wt = pro_df.head(1536)
    sig_d = defaultdict(lambda: defaultdict(float))

    Wt.index = sorted(list(snvs_signature_analyzer()))
    order_triplet = order_sigprofiler_snvs()
    sigs = Wt.columns
    for i, row in Wt.iterrows():
        for sig in sigs:
            colapsed = '{}[{}>{}]{}'.format(i[1], i[2], i[5], i[3])
            sig_d[sig][colapsed] += row[sig]

    collapsed_df = pd.DataFrame.from_dict(sig_d)
    collapsed_df = collapsed_df.loc[order_triplet]
    collapsed_df_normalized = collapsed_df / collapsed_df.sum()

    return collapsed_df_normalized


# this will check cosinus similarities with previously defined signatures and rename then accordingly
def assess_similar_signatures(d, exp_df, type_mut, good_sigs, outpath, ttype):
    # get PCAWG signatures to evaluate ML
    if type_mut == 'snvs':
        d_pcawg, index = pcawg_canonical_snvs()

    elif type_mut == 'dbs':
        d_pcawg, index = pcawg_canonical_dbs()
        index_f = [x.replace('>', '_') for x in index]
        index = index_f

    elif type_mut == 'indels':
        d_pcawg, index = pcawg_canonical_indels()

    # load processes found active in these samples
    d.columns = good_sigs

    original_sigs = d.columns.tolist()  # list of original signatures

    d = d / d.sum()

    fixed_cols = ['{}'.format(c.split('x')[1]) for c in d.columns]
    d.columns = fixed_cols
    d.index = index
    exp_df.index = original_sigs

    # find similar signatures to those reported in SigProfiler
    cos_sim = defaultdict(dict)
    for ix, col in enumerate(d.columns):
        vec1 = d[col].tolist()
        for s, vec2 in d_pcawg.items():
            c = 1 - round(spdist.cosine(vec1, vec2), 7)
            cos_sim[col][s] = c

    # select those with higher similarity per each of the signatures
    cos_df = pd.DataFrame(cos_sim)
    index_max = cos_df.idxmax()
    vals = cos_df.max()
    similar = pd.DataFrame(list(zip(index_max, vals)))

    similar.index = cos_df.columns

    # decide equivalent signatures with this 0.85 cutoff. Otherwise, give a new name
    new_sig = 1
    dict_equi = defaultdict(str)

    for i, row in similar.iterrows():
        if (row[1] > 0.85) & ('toremove' not in i):
            dict_equi[i] = '{}_{}_{}_{}'.format(new_sig, row[0], round(float(row[1]), 6), i.split('_')[1])
        else:
            if 'toremove' not in i:
                dict_equi[i] = '{}_{}'.format(new_sig, i.split('_')[1])
            else:
                dict_equi[i] = '{}_{}_{}_{}'.format(new_sig, row[0], new_sig, i)
        new_sig += 1

    cos_df_cols = cos_df.columns
    cos_df.columns = [dict_equi[c] for c in cos_df_cols]

    # plot cosinus similarity matrix
    fig, ax = plt.subplots(1, 1, figsize=(50, 20))
    config_params(7)
    sns.heatmap(cos_df.T, annot=True, fmt='g', cmap='YlGnBu')
    plt.savefig('{}/processes/{}/{}.heatmap.png'.format(outpath, ttype, ttype), dpi=300)
    plt.close()

    d.columns = [dict_equi[c] for c in cos_df_cols]

    return d


def process_exposures(exp_df, d, outpath, ttype):
    exp_df.index = d.columns
    # check samples just in case they are still in the full nomenclature
    cols = [c.split('_')[1] if 'post' in c else c for c in exp_df.columns]
    exp_df.columns = cols

    if ('Pan_' in ttype) or ('Pan.' in ttype) or (ttype == 'Pan'):
        cols = [c.split('.')[1] for c in exp_df.columns]
        exp_df.columns = cols

    exp_df.to_csv('{}/exposures/{}/{}.exposures.tsv'.format(outpath, ttype, ttype),
                  sep='\t', index=True, header=True)

    norm = exp_df / exp_df.sum()
    config_params(7)
    heatmap = sns.clustermap(norm.fillna(0), cmap='YlGnBu', figsize=(15, 10))
    heatmap.savefig('{}/exposures/{}/{}.heatmap.png'.format(outpath, ttype, ttype), dpi=300)
    plt.close()


# it will check for similar signatures, process processes and exposures
def signature_evaluation(d, exp_df, good_sigs, ttype, outpath, type_mut):
    # create the path if it is not found yet
    pathlib.Path('{}/processes/{}/'.format(outpath, ttype)).mkdir(parents=True, exist_ok=True)
    pathlib.Path('{}/exposures/{}/'.format(outpath, ttype)).mkdir(parents=True, exist_ok=True)

    outpath_process = '{}/processes/{}/{}.processes.tsv'.format(outpath, ttype, ttype)

    # check similar signatures
    d = assess_similar_signatures(d, exp_df, type_mut, good_sigs, outpath, ttype)
    d.to_csv(outpath_process, sep='\t', index=True, header=True)

    # Prepare the order for plotting
    if type_mut == 'snvs':
        d = d.loc[order_to_plot_snvs()]

    if type_mut == 'dbs':
        d = d.loc[order_to_plot_dbs()]

    if type_mut == 'indels':
        d = d.loc[order_to_plot_indel()]

    for s, sig in d.iteritems():
        if type_mut == 'snvs':
            plot_snvs(list(sig), s, outpath, ttype)

        elif type_mut == 'dbs':
            plot_dbs(list(sig), s, outpath, ttype)

        elif type_mut == 'indels':
            plot_indel(list(sig), s, outpath, ttype)

    # do heatmap and save exposures
    process_exposures(exp_df, d, outpath, ttype)


# extract info from SignatureAnalyzer run
def signature_matrix_siganalyzer(f, best_run):
    temp_dir = '{}/temp'.format(os.path.dirname(f))
    sample = os.path.basename(os.path.dirname(f))

    file_processes = "{}/{}.{}.processes.tsv".format(temp_dir, sample, best_run)

    df_processes = pd.read_csv(file_processes, sep='\t')
    df_processes.columns = [np.arange(1, len(df_processes.T) + 1)]

    sig_d = defaultdict(lambda: defaultdict(float))

    # merge channels if it belongs to the full composite
    if '_full' in f:

        df_processes.index = sorted(list(snvs_signature_analyzer()))
        order_triplet = order_sigprofiler_snvs()
        sigs = df_processes.columns
        for i, row in df_processes.iterrows():
            for sig in sigs:
                colapsed = '{}[{}>{}]{}'.format(i[1], i[2], i[5], i[3])
                sig_d[sig][colapsed] += row[sig]

        collapsed_df = pd.DataFrame.from_dict(sig_d)
        collapsed_df = collapsed_df.loc[order_triplet]

    else:
        collapsed_df = df_processes

    normalized_W = collapsed_df / collapsed_df.sum()
    normalized_W.to_csv('{}/{}.processes.tsv'.format(os.path.dirname(f), sample), sep='\t', index=False)
    collapsed_df.to_csv('{}/{}.processes_full.tsv'.format(os.path.dirname(f), sample), sep='\t', index=False)

    file_exposures = "{}/{}.{}.exposures.tsv".format(temp_dir, sample, best_run)
    df_exposures = pd.read_csv(file_exposures, sep='\t')
    df_exposures.to_csv('{}/{}.exposures.tsv'.format(os.path.dirname(f), sample), sep='\t', index=False)


# merge channels in bias analysis
def merge_channels_bias(d):
    dmerge = defaultdict(list)

    for c in d.columns:
        vec = d[c]
        mean_vec = []

        # get the mean of both contributions
        for lag, lead in zip(vec[0::2], vec[1::2]):
            mean_vec.append((lag + lead) / 2)

            dmerge[c] = mean_vec

    dmerged = pd.DataFrame(dmerge)
    dmerged_norm = dmerged / dmerged.sum()

    return dmerged_norm


# evaluate signatures when analysing asymmetries
def signature_evaluation_bias(d, exp_df, good_sigs, ttype, outpath, type_mut, bias):
    # create the path if it is not found yet
    pathlib.Path('{}/processes/{}/'.format(outpath, ttype)).mkdir(parents=True, exist_ok=True)
    pathlib.Path('{}/exposures/{}/'.format(outpath, ttype)).mkdir(parents=True, exist_ok=True)

    # labels of different biases
    d_bias = {'Replication': ("lagging", "leading"),
              'Transcription': ("transcribed", "untranscribed"),
              }

    labels_bias = d_bias[bias]

    dmerged_norm = merge_channels_bias(d)
    fixed_cols = ['{}'.format(c.split('x')[1]) for c in d.columns]

    dmerged_norm.columns = fixed_cols

    d_norm = assess_similar_signatures(dmerged_norm, exp_df, type_mut, good_sigs, outpath, ttype)

    d.columns = d_norm.columns

    # Plotting each of the process!!
    d_first_condition = d.iloc[::2]
    d_second_condition = d.iloc[1::2]
    d_second_condition.index = d_norm.index.tolist()
    d_first_condition.index = d_norm.index.tolist()

    toax = []
    for ord in d_norm.index.tolist():
        toax.append('{} {}'.format(ord, labels_bias[0]))
        toax.append('{} {}'.format(ord, labels_bias[1]))

    d.index = toax
    dnorm = d / d.sum()

    outpath_process = '{}/processes/{}/{}.processes.tsv'.format(outpath, ttype, ttype)
    dnorm.to_csv(outpath_process, sep='\t', index=True, header=True)

    if type_mut == 'snvs':
        d_clustered = d_first_condition.loc[order_to_plot_snvs()]
        d_unclustered = d_second_condition.loc[order_to_plot_snvs()]
        plot_bias_snv(d_clustered, d_unclustered, outpath, ttype, labels_bias[0], labels_bias[1])

    if type_mut == 'dbs':
        d_clustered = d_first_condition.loc[order_to_plot_dbs()]
        d_unclustered = d_second_condition.loc[order_to_plot_dbs()]
        plot_bias_dbs(d_clustered, d_unclustered, outpath, ttype, labels_bias[0], labels_bias[1])

    if type_mut == 'indels':
        d_clustered = d_first_condition.loc[order_to_plot_indel()]
        d_unclustered = d_second_condition.loc[order_to_plot_indel()]
        plot_bias_indel(d_clustered, d_unclustered, outpath, ttype, labels_bias[0], labels_bias[1])

    process_exposures(exp_df, d, outpath, ttype)


# get stabilities and merge all files after the fitting
def create_summary_sigprofiler(original_genome, results_path):

    ttype_path = os.path.basename(original_genome).split('.dlm')[0]
    outpath = '{}/{}/'.format(results_path, ttype_path)

    summary_file = '{}/SigProfiler_summary.pckl.gz'.format(outpath)

    if not os.path.isfile(summary_file):

        stability_processes = 'processesStabAvg'
        processes_file = 'processes'
        exposures_file = 'exposures_fitting'
        avg_reconstruction = 'avgReconstructionError'

        dic_stats = defaultdict(lambda: defaultdict(dict))
        sigs_examined = []

        for nsigs_stabavg in glob('{}/processesStabAvg_*'.format(outpath)):

            sig_n = int(os.path.basename(nsigs_stabavg).split('_')[1])
            sigs_examined.append(sig_n)

            file_stability_processes = os.path.join(outpath, stability_processes + "_" + str(sig_n))
            file_processes = os.path.join(outpath, processes_file + "_" + str(sig_n))
            file_exposures = os.path.join(outpath, exposures_file + "_" + str(sig_n))
            file_avgreconserror = os.path.join(outpath, avg_reconstruction + "_" + str(sig_n))

            dic_stats[sig_n]['ProcessFile'] = file_processes
            dic_stats[sig_n]['ExposuresFile'] = file_exposures

            # read the stability of the signatures
            df_signature_stab = pd.read_csv(file_stability_processes, sep='\t')

            # add all the values to a dictionary
            dic_stats[sig_n]['stabilitySigs'] = df_signature_stab.values[0]

            # get median and mean stability
            median_sig = df_signature_stab.median(axis=1)
            mean_sig = df_signature_stab.mean(axis=1)

            # keep it in dictionary
            dic_stats[sig_n]['meanStab'] = mean_sig
            dic_stats[sig_n]['medianStab'] = median_sig
            dic_stats[sig_n]['nameSigs'] = df_signature_stab.columns.tolist()

            # read the stability of the signatures
            df_avgreconserror = pd.read_csv(file_avgreconserror, sep='\t', header=None)

            dic_stats[sig_n]['FroError'] = float(df_avgreconserror.loc[0])

        pickle_out = '{}/SigProfiler_summary.pckl.gz'.format(outpath)
        pickle.dump(dict(dic_stats), gzip.open(pickle_out, 'wb'))


# plot to decide the number of active signatures according to SigProfiler
def extract_best_signatures_sigprofiler(original_genome, results_path, final_path, similarity_cutoff):

    os.makedirs(final_path, exist_ok=True)

    ttype = os.path.basename(original_genome).split('.')[0]
    ttype_path = os.path.basename(original_genome).split('.dlm')[0]
    current_path = '{}/{}/'.format(results_path, ttype_path)

    d = pickle.load(gzip.open('{}/SigProfiler_summary.pckl.gz'.format(current_path)))

    sigs_examined = list(d.keys())

    fig, ax = plt.subplots(1, 1, figsize=(16, 3))

    min_sig = min(sigs_examined)
    max_sig = max(sigs_examined)

    plt.title(ttype)

    plt.xticks([i for i in range(min_sig, max_sig)])
    plt.xlabel('Signatures')

    ax2 = ax.twinx()
    ax3 = ax.twinx()

    ax3.spines["right"].set_position(("axes", 1.06))

    ax.hlines(0.95, min_sig, max_sig, linestyles='--', color='grey', alpha=0.5)
    ax.hlines(0.90, min_sig, max_sig, linestyles='--', color='grey', alpha=0.5)
    ax.hlines(0.85, min_sig, max_sig, linestyles='--', color='grey', alpha=0.5)

    for ix, (sig, dic) in enumerate(d.items()):

        # stabilities
        ax.scatter(sig, dic['medianStab'], c='darkred', label='Median')
        ax.scatter(sig, dic['meanStab'], c='orange', label='Mean')
        ax2.scatter(sig, dic['FroError'], c='darkblue', label='FroErr')
        nhigher = len([i for i in dic['stabilitySigs'] if i > similarity_cutoff])
        ax3.scatter(sig, nhigher, c='purple', label='HighRep Sigs')

        if ix == 0:
            ax.legend(bbox_to_anchor=(-0.05, 1.10))
            ax2.legend(bbox_to_anchor=(-0.05, 0.9))
            ax3.legend(bbox_to_anchor=(-0.05, 0.7))

    ax.set_ylim(0.50, 1)
    ax.set_ylabel('Stability')
    ax2.set_ylabel('Reconstruction error')
    ax3.set_ylabel('Number of signatures')

    plt.savefig('{}/{}.png'.format(final_path, ttype_path), dpi=300, bbox_inches="tight")

    plt.close()


# main function for SigProfiler
def postprocessing_sigprofiler(dic_signatures_active, results_path, similarity_cutoff, outpath):

    for ttype, nsigs_active in dic_signatures_active.items():

        type_mut = ttype.split('.')[-1]

        print("Doing {} with {} sigs active, {}".format(ttype, nsigs_active, type_mut))

        current_path = '{}/{}/'.format(results_path, ttype)
        d = pickle.load(gzip.open('{}/SigProfiler_summary.pckl.gz'.format(current_path)))

        processed_file = d[nsigs_active]['ProcessFile']
        exposures_file = d[nsigs_active]['ExposuresFile']
        name_sigs = d[nsigs_active]['nameSigs']

        pro_df = pd.read_csv(processed_file, sep='\t')
        exp_df = pd.read_csv(exposures_file, sep='\t')

        # label those unstable
        good_sigs = ['{}_{}'.format(i, round(v, 2)) if v > similarity_cutoff else
                     '{}_toremove_{}_'.format(i, round(v, 2))
                     for i, v in zip(name_sigs, d[nsigs_active]['stabilitySigs'])]

        outpath_final = '{}/{}/'.format(outpath, type_mut)

        signature_evaluation(pro_df, exp_df, good_sigs, ttype, outpath_final, type_mut)


# main function for SignatureAnalyzer
def postprocessing_signatureanalyzer(summaries, outpath_init, bias):
    all_files = glob(summaries)
    for exposures_file in all_files:

        if bias != 'None':
            processed_file = exposures_file.replace('exposures', 'processes_full')
        else:
            processed_file = exposures_file.replace('exposures', 'processes')

        ttype = os.path.basename(exposures_file).split('.')[0]

        type_mut = os.path.basename(exposures_file).split('.')[1]
        outpath = '{}/{}/'.format(outpath_init, type_mut)

        if (os.path.isfile(processed_file)) & (os.path.isfile(exposures_file)):

            pro_df = pd.read_csv(processed_file, sep='\t')
            exp_df = pd.read_csv(exposures_file, sep='\t')
            exp_df.reset_index(drop=True, inplace=True)

            new_cols = ['x{}_1'.format(c) for c in pro_df.columns]
            pro_df.columns = new_cols
            good_sigs = new_cols

            # collapse if we do the full approach
            if '_full' in ttype:
                pro_df = return_1536_to_96_channels(pro_df)

            if bias != "None":
                signature_evaluation_bias(pro_df, exp_df, good_sigs, ttype, outpath, type_mut, bias)
            else:
                signature_evaluation(pro_df, exp_df, good_sigs, ttype, outpath, type_mut)
        else:
            print('processes and exposures not found')
