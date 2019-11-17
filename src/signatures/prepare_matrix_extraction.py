# Import modules
import os
import sys
from collections import defaultdict

import click
import pandas as pd
import pathlib
from pybedtools import BedTool

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
scripts_path = os.path.abspath(os.path.join(__file__,  '..', '..'))
sys.path.insert(0, scripts_path)

from sign_utils.metadata import return_metadata
from sign_utils.orders import snvs_signature_analyzer, order_muts
from config import Conf


# Global variables
rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


# put mutations in the C,T framework
def fix_extended_context(df):

    if df['REF'] not in 'CT':
        ext_original = df['EXT_MUT']
        ext = '{}{}{}{}{}{}'.format(
            rev[ext_original[-2]],
            rev[ext_original[-3]],
            rev[ext_original[-4]],
            rev[ext_original[-5]],
            rev[ext_original[-6]],
            rev[ext_original[-1]]
        )
        toreturn = ext
    else:
        toreturn = df['EXT_MUT']

    return toreturn


def return_ordered_matrix(df, order, file, col='VARIANT_CLASS'):

    samples_dict = defaultdict(dict)
    d_metadata, dic_secondary_fixed = return_metadata()

    for sample, data in df.groupby(by='SAMPLE'):

        dic_count = data[col].value_counts().to_dict()
        for i in order:
            if ('Pan.' in file) or ('PanHypermutated' in file) or ('PanNotHypermutated' in file):
                samples_dict['{}.{}'.format(d_metadata[sample], sample)][i] = dic_count.get(i, 0)
            else:
                samples_dict[sample][i] = dic_count.get(i, 0)
    matrix = pd.DataFrame.from_dict(samples_dict)

    return matrix


# prepare the matrix for SignatureAnalyzer. If we are looking at Pan, it will include the ttype info
def prepareMatrixExtraction(t):

    snv_file, outpath = t
    # get also the indels and dbs files
    namef = os.path.basename(snv_file).split('.snvs.gz')[0]
    indel_file = snv_file.replace('snvs', 'indels')
    dbs_file = snv_file.replace('snvs', 'dbs')

    df = pd.read_csv(snv_file, sep='\t', low_memory=True, keep_default_na=False)

    pathlib.Path('{}'.format(outpath)).mkdir(parents=True, exist_ok=True)

    if len(df) > 0:

        # SNVs PENTAMER
        df['EXT_MUT'] = df.apply(lambda x: '{}{}'.format(x['EXTENDED'], x['ALT']), axis=1)
        df['FIX_EXT'] = df.apply(fix_extended_context, axis=1)

        order = snvs_signature_analyzer()
        matrix_pentamers = return_ordered_matrix(df, order, snv_file, col='FIX_EXT')

        # SNVs TRIPLET
        order = order_muts('snv')

        matrix_snvs_triplet = return_ordered_matrix(df, order, snv_file)

        # Indels
        matrix_indels = []
        if os.path.isfile(indel_file):
            df = pd.read_csv(indel_file, sep='\t', low_memory=True, keep_default_na=False)
            order = order_muts('indels')
            matrix_indels = return_ordered_matrix(df, order, snv_file)

        # DBS
        matrix_dbs = []
        if os.path.isfile(dbs_file):
            df = pd.read_csv(dbs_file, sep='\t', low_memory=True, keep_default_na=False)
            order = order_muts('dbs')
            matrix_dbs = return_ordered_matrix(df, order, snv_file)

        # FULL for SignatureAnalyzer
        total = [matrix_pentamers, matrix_indels, matrix_dbs]
        outsnv = pd.concat(total)

        # save outputs
        outsnv.fillna(0).to_csv('{}/{}_full.snvs.dlm'.format(outpath, namef), sep='\t', index=False)
        matrix_dbs.fillna(0).to_csv('{}/{}.dbs.dlm'.format(outpath, namef), sep='\t', index=False)
        matrix_indels.fillna(0).to_csv('{}/{}.indels.dlm'.format(outpath, namef), sep='\t', index=False)
        matrix_snvs_triplet.fillna(0).to_csv('{}/{}.snvs.dlm'.format(outpath, namef), sep='\t', index=False)


# classify mutations in replication counts
def classifyReplication(df):

    # reference
    ref = 'CT'

    # this means canonical is the leading
    if df['is_left'] == "1":
        if df['REF'] in ref:
            var_class = '{}_leading'.format(df['VARIANT_CLASS'])
        else:
            var_class = '{}_lagging'.format(df['VARIANT_CLASS'])

    # this means canonical is the lagging
    elif df['is_right'] == "1":
        if df['REF'] in ref:
            var_class = '{}_lagging'.format(df['VARIANT_CLASS'])
        else:
            var_class = '{}_leading'.format(df['VARIANT_CLASS'])
    else:
        var_class = 'unclassified'

    return var_class


def prepareMatrixExtraction_replication(tup):

    f, outpath = tup
    namef = os.path.basename(f).split('.snvs.gz')[0]
    d_metadata, dic_secondary_fixed = return_metadata()

    for type_var in ['snvs', 'dbs', 'indels']:

        var_file = f.replace('snvs', type_var)
        # regions obtained from Haradvala
        timing_file = Conf['replication_file']

        df_time = pd.read_csv(timing_file, sep='\t')

        # get only regions with high enough slope, similar to Haradhvala et al
        RT_regions = df_time[df_time['rt'] > 250].sort_values(by=['chr', 'st'])

        df = pd.read_csv(var_file, sep='\t')
        df['POS-1'] = df['POS'] - 1

        # intersect with RT bins
        bedmut = BedTool.from_dataframe(
            df[['CHROM', 'POS-1', 'POS', 'REF', 'ALT', 'TRIPLET', 'SAMPLE', 'VARIANT_CLASS']]
        )
        intersected = bedmut.intersect(BedTool.from_dataframe(RT_regions), wao=True, sorted=True).to_dataframe(
            names=[
                'CHROM', 'POS-1', 'POS', 'REF', 'ALT', 'TRIPLET', 'SAMPLE', 'VARIANT_CLASS',
                'chr', 'st', 'en', 'is_left', 'is_right', 'txplus', 'txminus', 'rt', 'overlapp'
            ], low_memory=False
        )
        # classify variants according to leading and lagging
        intersected['class_var'] = intersected.apply(classifyReplication, axis=1)

        # remove unclassified regions
        intersected = intersected[intersected['class_var'] != 'unclassified']

        if type_var == 'snvs':
            to_order = 'snv'
        else:
            to_order = type_var

        order = order_muts(to_order)
        samples_dict = defaultdict(dict)

        for sample, data in intersected.groupby(by='SAMPLE'):

            if ('Pan.' in var_file) or ('PanHypermutated' in var_file) or ('PanNotHypermutated' in var_file):
                sample = '{}.{}'.format(d_metadata[sample], sample)

            dic_count = data['class_var'].value_counts().to_dict()

            for i in order:
                order_lag = '{}_lagging'.format(i)
                order_lead = '{}_leading'.format(i)
                samples_dict[sample][order_lag] = dic_count.get(order_lag, 0)
                samples_dict[sample][order_lead] = dic_count.get(order_lead, 0)

        matrix = pd.DataFrame.from_dict(samples_dict)
        out = '{}/{}_replication.{}.dlm'.format(outpath, namef, type_var)

        matrix.to_csv(out, sep='\t', header=True, index=False)


# classify mutations in transcribed counts
def classifyTranscibed(df):

    # reference
    ref = 'CT'

    # this means canonical is the leading
    if df['strand'] == "+":
        if df['REF'] in ref:
            var_class = '{}_transcribed'.format(df['VARIANT_CLASS'])
        else:
            var_class = '{}_untranscribed'.format(df['VARIANT_CLASS'])

    # this means canonical is the lagging
    elif df['strand'] == "-":
        if df['REF'] in ref:
            var_class = '{}_untranscribed'.format(df['VARIANT_CLASS'])
        else:
            var_class = '{}_transcribed'.format(df['VARIANT_CLASS'])
    else:
        var_class = 'unclassified'

    return var_class


def prepareMatrixExtraction_transcription(tup):

    f, outpath = tup

    namef = os.path.basename(f).split('.snvs.gz')[0]
    d_metadata, dic_secondary_fixed = return_metadata()

    for type_var in ['snvs', 'dbs', 'indels']:

        var_file = f.replace('snvs', type_var)

        # regions obtained from Haradvala
        transcription_file = Conf['transcription_file']

        df_transcribed = pd.read_csv(transcription_file, sep='\t', names=['chr', 'st', 'en', 'strand'])
        df_transcribed.sort_values(by=['chr', 'st'], inplace=True)

        df = pd.read_csv(var_file, sep='\t')

        df['POS-1'] = df['POS'] - 1

        # intersect with RT bins
        bedmut = BedTool.from_dataframe(
            df[['CHROM', 'POS-1', 'POS', 'REF', 'ALT', 'TRIPLET', 'SAMPLE', 'VARIANT_CLASS']]
        )
        intersected = bedmut.intersect(BedTool.from_dataframe(df_transcribed), wao=True, sorted=True).to_dataframe(
            names=[
                'CHROM', 'POS-1', 'POS', 'REF', 'ALT', 'TRIPLET', 'SAMPLE',
                'VARIANT_CLASS', 'chr', 'st', 'en', 'strand', 'overlapp'
            ], low_memory=False
        )

        # classify variants according to transcribed or untranscribed
        if len(intersected) > 0:

            intersected['class_var'] = intersected.apply(classifyTranscibed, axis=1)

            # remove unclassified regions
            intersected = intersected[intersected['class_var'] != 'unclassified']

            if type_var == 'snvs':
                to_order = 'snv'
            else:
                to_order = type_var

            order = order_muts(to_order)
            samples_dict = defaultdict(dict)

            for sample, data in intersected.groupby(by='SAMPLE'):

                if ('Pan.' in var_file) or ('PanHypermutated' in var_file) or ('PanNotHypermutated' in var_file):
                    sample = '{}.{}'.format(d_metadata[sample], sample)

                dic_count = data['class_var'].value_counts().to_dict()
                for i in order:
                    order_lag = '{}_transcribed'.format(i)
                    order_lead = '{}_untranscribed'.format(i)
                    samples_dict[sample][order_lag] = dic_count.get(order_lag, 0)
                    samples_dict[sample][order_lead] = dic_count.get(order_lead, 0)

            matrix = pd.DataFrame.from_dict(samples_dict)
            out = '{}/{}_transcription.{}.dlm'.format(outpath, namef, type_var)
            matrix.to_csv(out, sep='\t', header=True, index=False)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    '-c', '--merged-cohort',
    required=True,
    help='Merged cohort (for instance data/hartwig/origin/Adrenal.snvs.gz)',
    type=click.Path()
)
@click.option(
    '-o', '--output-dir',
    required=True,
    help='Output directory with the results (for instance data/hartwig/signatures/matrix/)',
    type=click.Path()
)
def matrix_creation(merged_cohort, output_dir):
    tup = (merged_cohort, output_dir)
    # prepare canonical matrix
    prepareMatrixExtraction(tup)

    # prepare replication/transcription matrix
    prepareMatrixExtraction_replication(tup)
    prepareMatrixExtraction_transcription(tup)


if __name__ == '__main__':
    matrix_creation()
