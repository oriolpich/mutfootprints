# Import modules
import os
import sys
from collections import defaultdict

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.neighbors.kde import KernelDensity
from scipy.signal import find_peaks


# Global variables
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


# ===
# This will generate a file with the clonal structure (it does not have to be perfect given than MT will already
# perform a E/MX approach )
# ===


def new_ccf(df, purity):

    vaf = df['VAR_COUNTS'] / (df['REF_COUNTS'] + df['VAR_COUNTS'])

    # definition of CCF
    ccf = vaf * (purity * df['TOTAL_CN'] + (1 - purity) * df['NORMAL_CN']) / purity
    return ccf


def do_clustering(outdf, pur, best_band, lower_cuttoff, upper_cuttoff):

    converged = 0

    # get the ccf value with the purity correction
    outdf['new_ccf'] = outdf.apply(new_ccf, args=(pur,), axis=1)

    # remove extreme cases
    ccf_list = [x for x in outdf[outdf['new_ccf'] < 2.8]['new_ccf'].tolist() if x > 0]
    max_ccf = np.amax(outdf['new_ccf'])

    if max_ccf < 2.8:
        upbound = max_ccf
    else:
        upbound = 2.8

    # do the log2 of each of the ccf values
    ccf = [np.log2(x) for x in ccf_list if x > 0]

    X = np.array(ccf).reshape(-1, 1)

    kde = KernelDensity(kernel='gaussian', bandwidth=best_band).fit(X)

    grid2 = np.linspace(np.amin(ccf_list), upbound, num=150).reshape(-1, 1)
    grid2 = np.array([np.log2(x) for x in grid2])
    flat_array = grid2.flatten()

    log_density = kde.score_samples(grid2)
    density = np.exp(log_density)

    # find the maximum peaks
    number_components = len(find_peaks(density, height=0.1)[0])
    gmm = GaussianMixture(n_components=number_components, max_iter=2000).fit(X)

    cluster_assign = defaultdict(list)

    for ix, prob in enumerate(np.argmax(gmm.predict_proba(X), axis=1)):
        cluster_assign[prob].append(X[ix])

    good_center_ccf = []
    for clust, vals in cluster_assign.items():

        center = np.median(vals)
        if lower_cuttoff <= center < upper_cuttoff:
            converged = pur
            good_center_ccf = center

    return converged,  good_center_ccf, cluster_assign, flat_array, density


# plot the ccf distribution and clusters in normal scale
def plot_ccf(cluster_assign, flat_array, density, outfile_plot):

    fig, ax = plt.subplots(1, 1, figsize=(2, 2))
    for clust, vals in cluster_assign.items():
        center = np.median(vals)
        ax.text(center + 0.1, 1, round(2 ** center, 3))
        ax.vlines(center, 0, 1, color='red')

    ax.plot(flat_array, density, c='#31a354', lw=2)
    plt.xlabel('ccf')
    xtickslocs = ax.get_xticks()
    ax.set_xticklabels([2 ** i for i in xtickslocs])
    plt.tight_layout()
    plt.savefig(outfile_plot)
    plt.close()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--mutation-file', required=True, help='Input mutation file', type=click.Path())
@click.option('-o', '--outpath', required=True, help='Output directory with the results', type=click.Path())
def clustering_ccf(mutation_file, outpath):

    os.makedirs(outpath, exist_ok=True)

    # define outfiles
    outfile = '{}{}.subclonal.txt'.format(outpath, os.path.basename(mutation_file).split('.')[0])
    outfile_plot = '{}{}.subclonal.png'.format(outpath, os.path.basename(mutation_file).split('.')[0])

    df = pd.read_csv(mutation_file, sep='\t')

    # select only those cases where REF_COUNTS > 0
    outdf = df[df['REF_COUNTS'] > 0]

    # get purity
    purity = outdf.iloc[0]['PURITY']

    # Round purity
    totrack = [purity]

    # explore purity range
    for i in np.arange(0.01, 0.06, 0.01):
        totrack.append(purity + i)
        totrack.append(purity - i)

    # cuttofs to say convergence
    upper_cuttoff = np.log2(1.07)
    lower_cuttoff = np.log2(0.93)

    # new purity converge
    converged = 0

    # hardcoded!
    best_band = 0.09

    # do the clustering and check if we converge with a cluster at ccf-1
    for ix, pur in enumerate(totrack):

        if converged == 0:
            converged, good_center_ccf, cluster_assign, flat_array, density = do_clustering(
                outdf, pur, best_band, lower_cuttoff, upper_cuttoff
            )

    if converged > 0:

        plot_ccf(cluster_assign, flat_array, density, outfile_plot)

        # get subclonal status
        subclonals = 0
        keep = []
        for center, vals in cluster_assign.items():
            if np.median(vals) < good_center_ccf:
                keep.append(((2 ** np.median(vals)) * converged, len(vals)))
                subclonals += len(vals)

        # get outvalues for MutationalTime
        nonsubclonal = len(df) - subclonals
        with open(outfile, 'wt') as outfile_name:
            out = 'proportion\tn_ssms\n'
            outfile_name.write(out)
            out = '{}\t{}\n'.format(converged, nonsubclonal)
            outfile_name.write(out)
            for subs in keep:
                conv, vals = subs
                out = '{}\t{}\n'.format(conv, vals)
                outfile_name.write(out)
    else:
        # put them on the blacklist
        outpath = '{}/blacklist/'.format(outpath)
        os.makedirs(outpath, exist_ok=True)
        outfile_plot = '{}{}.subclonal.png'.format(outpath, os.path.basename(f).split('.')[0])

        plot_ccf(cluster_assign, flat_array, density, outfile_plot)


if __name__ == '__main__':
    clustering_ccf()
