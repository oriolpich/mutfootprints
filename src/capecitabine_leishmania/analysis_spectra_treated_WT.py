# Import modules
import os
import sys
from glob import glob
from itertools import islice
from collections import defaultdict, Counter
import io

import click
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

os.environ['BGDATA_BUILDS'] = "src/utils/bgdata_bgreference.txt"
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]

from bgreference import refseq

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..',))
sys.path.insert(0, scripts_path)

print(scripts_path)
from utils.orders import order_muts
from utils.signature_plots import config_params, chunks


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def slicing_window(seq, n=3):
    it = iter(seq)
    result = ''.join(islice(it, n))

    if len(result) == n:
        yield result

    for elem in it:
        result = result[1:] + elem
        yield result


def composition_WG():

    all_chroms = glob('data/capecitabine_leishmania/bgreference/leish/*.txt')
    dic_count = Counter()
    for f in tqdm(all_chroms):
        if 'CACT' not in f:
            with open(f, 'rt') as infile:
                for line in infile:
                    line = line.rstrip()
                    sliced = list(slicing_window(line))
                    dic_count = dic_count + Counter(sliced)

    return dic_count


composition_genome = composition_WG()


def create_snv_class(df):

    pyr = ['C', 'T']
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    x = df['TRIPLET']

    if x[1] in pyr:
        out = '{}[{}>{}]{}'.format(x[0], x[1], df['ALT'], x[2])
    else:
        out = '{}[{}>{}]{}'.format(rev[x[2]], rev[x[1]], rev[df['ALT']], rev[x[0]])
    return out


# vcf parser
def vcf_reader(path):
    with open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(str.join(os.linesep, lines)),
        dtype={
            '#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str
        }, sep='\t', low_memory=False,
    ).rename(columns={'#CHROM': 'CHROM'})


def plot_single_sig(sig, order_plot, outfile):

    config_params(3)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(3.2, 1),
                            gridspec_kw={'height_ratios': [1, 9]})

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

    axs[0].set_title(outfile)

    x = [i for i in range(len(vals))]

    axs[1].axhline(y=0.05, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.1, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)
    axs[1].axhline(y=0.15, xmin=-1, xmax=96, lw=0.6, color='grey', alpha=0.2)

    axs[1].bar(x, vals, color=colors, width=0.8, linewidth=0, align='center')
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(['{}{}{}'.format(a[0], a[2], a[-1]) for a in order_plot],
                           verticalalignment="center", ha='center', rotation=90, fontsize=2,
                           color='grey')

    plt.tight_layout()
    plt.xlim(-1, 96)

    # axs.spines['right'].set_visible(False)
    axs[1].spines['top'].set_visible(False)
    axs[1].set_ylabel('Relative Probability')

    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=0.5)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)
    plt.savefig('figures/FU5_{}.svg'.format(outfile))


def do_plot(dic_variants, outfile):

    dic_norm = defaultdict(float)
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    for mut, count in dic_variants.items():
        triplet = '{}{}{}'.format(mut[0], mut[2], mut[-1])
        rev_triplet = '{}{}{}'.format(rev[mut[-1]], rev[mut[2]], rev[mut[0]])
        total_count = composition_genome[triplet] + composition_genome[rev_triplet]
        dic_norm[mut] = count / total_count

    order_plot = order_muts('snv')
    sigs = pd.DataFrame(dic_norm.items())
    sigs = sigs.set_index(0)
    sigs = sigs.loc[order_plot].fillna(0)
    sigs = sigs / sigs.sum()

    print(sigs[1].tolist())

    plot_single_sig(sigs[1].tolist(), order_plot, outfile)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--vcf-file', required=True, help='Input VCF file', type=click.Path())
def analyze_experiment(vcf_file):

    df = vcf_reader(vcf_file)
    df['POS-1'] = df['POS']-1

    # remove non canonical chromosomes
    df['CAN'] = df['CHROM'].apply(lambda x: 'RMV' if 'CACT' in x else 'PASS')
    df = df[df['CAN'] == 'PASS']

    df['TRIPLET'] = df.apply(lambda x: refseq('leish', x['CHROM'], x['POS-1'], 3, release=None), axis=1)

    # select whether we have SNVs or others
    df['len_alt'] = df['ALT'].str.len()

    # number of characters in ref
    df['len_ref'] = df['REF'].str.len()

    # first classification between SNV and others
    df['TYPE'] = df.apply(
        lambda x: 'SNV' if (x['len_alt'] == 1) & (x['len_ref'] == 1) and
                           (x['ALT'] != '-') and (x['REF'] != '-') else 'INDEL', axis=1
    )

    df = df[df['TYPE'] == 'SNV']
    df['VARIANT_CLASS'] = df.apply(create_snv_class, axis=1)

    # get whether the mutation has happened in the WT or treated
    df['FILTER_WT'] = df['ERR174230.bam'].apply(lambda x: 'RES' if '0/0' in x else 'WT')

    # select only variants in the treated
    dic_variants = df[df['FILTER_WT'] == 'RES']['VARIANT_CLASS'].value_counts().to_dict()

    do_plot(dic_variants, 'experiment')

    # select only variants in the WT
    dic_variants = df[df['FILTER_WT'] == 'WT']['VARIANT_CLASS'].value_counts().to_dict()
    do_plot(dic_variants, 'WT')


if __name__ == '__main__':
    analyze_experiment()
