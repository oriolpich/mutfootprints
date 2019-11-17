# Import modules
import gzip
import os
import pickle
import sys

import click
import pandas as pd
from pybedtools import BedTool


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def load_CNAS():

    d = pickle.load(gzip.open('data/hartwig/clonal_structure/dic_cnas.pckl.gz'))
    return d


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--mutation-file', required=True, help='Input mutation file', type=click.Path())
@click.option('-o', '--outpath', required=True, help='Output directory with the results', type=click.Path())
def get_CNAS_ccf(mutation_file, outpath):

    dic_cnas = load_CNAS()

    # only the name of the tumoral sample
    name = os.path.basename(mutation_file).split('_')[1].split('.')[0]

    # full name
    fullname = os.path.basename(mutation_file).split('.')[0]
    outpath = '{}/{}.cna.gz'.format(outpath, fullname)

    cna_file = dic_cnas[name]
    purity_file = dic_cnas[name].replace('.cnv', '.purity')

    # read CNAs
    df = pd.read_csv(cna_file, sep='\t')
    df['#chromosome'] = df['#chromosome'].apply(lambda x: 'chr{}'.format(x))
    df['major_cn'] = round(df['baf'] * df['copyNumber'])
    df['minor_cn'] = round((1 - df['baf']) * df['copyNumber'])
    df.replace(-0.0, 0, inplace=True)
    filtered_df = df[(df['major_cn'] >= 0) & (df['minor_cn'] >= 0)]

    # read purity and gender
    pur = pd.read_csv(purity_file, sep='\t')
    purity = pur['#Purity'].tolist()[0]
    gender = pur['Gender'].tolist()[0].lower()

    # read the mutations
    df_muts = pd.read_csv(mutation_file, sep='\t')

    df_muts['POS-1'] = df_muts['POS'] - 1
    mut_bed = BedTool.from_dataframe(df_muts[['CHROM', 'POS-1', 'POS']])
    mut_cnas = BedTool.from_dataframe(
        filtered_df[['#chromosome', 'start', 'end', 'major_cn', 'minor_cn']]
    ).coverage(mut_bed, counts=True)
    dfc = mut_cnas.to_dataframe(names=[
        'chromosome', 'start', 'end', 'major_cn', 'minor_cn', 'n.snv_mnv'
    ])

    dfc['clonal_frequency'] = purity
    dfc['gender'] = gender

    dfc[['chromosome', 'start', 'end', 'major_cn', 'minor_cn', 'clonal_frequency', 'gender', 'n.snv_mnv']].to_csv(
        outpath, header=True, index=False, sep='\t', compression='gzip'
    )


if __name__ == '__main__':
    get_CNAS_ccf()
