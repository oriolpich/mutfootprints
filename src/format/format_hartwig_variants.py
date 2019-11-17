import io
import os
import sys
import gzip

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

import pandas as pd
from pybedtools import BedTool
import pybedtools
import numpy as np
from bgreference import hg19
from utils.indels_class import IndelsClassifier
from config import Conf


# vcf parser
def vcf_reader(path):
    with gzip.open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(str.join(os.linesep, lines)),
        dtype={
            '#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str
        }, sep='\t', low_memory=False,
    ).rename(columns={'#CHROM': 'CHROM'})


# get reads info from VCF file
def get_reads(df, last_column):
    list_format = df['FORMAT'].split(':')
    d = {list_format[i]: i for i in range(len(list_format))}

    colvaf_list = df[last_column].split(':')

    reference = int(colvaf_list[d['AD']].split(',')[0])
    variant = int(colvaf_list[d['AD']].split(',')[1])

    total_reads = reference + variant
    df['total_reads'] = total_reads
    df['ref_reads'] = reference
    df['var_reads'] = variant

    if total_reads > 0:
        df['VAF'] = df['var_reads'] / df['total_reads']
    else:
        df['VAF'] = np.nan

    return df


# create order format similar to SigProfiler
def create_snv_class(df):
    pyr = ['C', 'T']
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    x = df['TRIPLET']
    if x[1] in pyr:
        out = '{}[{}>{}]{}'.format(x[0], x[1], df['ALT'], x[2])
    else:
        out = '{}[{}>{}]{}'.format(rev[x[2]], rev[x[1]], rev[df['ALT']], rev[x[0]])
    return out


# Get the major CN in the tumor. If it has no change, then we put the normal cn
def get_major_cn(df):
    if df['overlapp'] == 1:
        major = round(np.float(df['MAJOR_CN_TEMP']))

    else:
        major = round(np.float(df['NORMAL_CN']))
    return major


# get mutations only in mappable regions
def get_mappable_regions(bed):
    path_mappable = Conf["mappability_file"]
    mapp_bed = BedTool(path_mappable)
    mappable_mut = bed.intersect(mapp_bed, wa=True)

    return mappable_mut


def load_files(mutation_file, cnvs_file, purity_file):

    # load cnvs file from purple output
    if os.path.isfile(cnvs_file) is not True:
        print('CNS file {} does not exist, exiting...'.format(cnvs_file))
        sys.exit()

    # read CNS file and store it in BedTool format
    df_cns = pd.read_csv(cnvs_file, sep='\t')
    cnv_bed = BedTool.from_dataframe(df_cns[['#chromosome', 'start', 'end', 'copyNumber', 'baf']])

    # load purity and gender
    if os.path.isfile(purity_file) is not True:
        print('Purity file file {} does not exist'.format(purity_file))
        sys.exit()

    df_purity = pd.read_csv(purity_file, sep='\t')
    purity_score = np.float(df_purity['#Purity'].tolist()[0])
    gender = df_purity['Gender'].tolist()[0]

    # read vcf
    df = vcf_reader(mutation_file)

    # get only canonical chromosomes
    wantedchroms = [str(i) for i in range(1, 23)]
    wantedchroms.append('Y')
    wantedchroms.append('X')

    # select only variants with the PASS filter and in the chromosomes we are interested in
    df = df[df['CHROM'].isin(wantedchroms)]
    df = df[df['FILTER'] == 'PASS']

    return df, cnv_bed, purity_score, gender


def indels_classification(indel_df):

    # Remove those cases where we have "," in the ALT
    indel_df = indel_df[~indel_df['ALT'].str.contains(',')]

    # set up first classification of variants
    indel_df['CLASS'] = indel_df.apply(
        lambda x: 'INS' if (x['len_alt'] > 1) & (x['len_ref'] == 1) else (
            'DEL' if (x['len_ref'] > 1) & (x['len_alt'] == 1) else (
                'DBS' if (x['len_alt'] == 2) & (x['len_ref'] == 2) else (
                    'MNV' if (x['len_alt'] == x['len_ref']) else (
                        'COMPLEX_INDELS')))), axis=1
    )

    complex_indels_df = indel_df[indel_df['CLASS'] == 'COMPLEX_INDELS']
    complex_indels_df['VARIANT_CLASS'] = 'COMPLEX_INDELS'
    mnv_indels_df = indel_df[indel_df['CLASS'] == 'MNV']
    mnv_indels_df['VARIANT_CLASS'] = 'MNV'

    # we won't be subclassifying these types
    toremove_class = ['MNV', 'COMPLEX_INDELS']

    # remove unwanted others
    indel_df = indel_df[~indel_df['CLASS'].isin(toremove_class)]

    # assing each class to the indels and dbs using PCAWG classification
    indel_df['VARIANT_CLASS'] = indel_df.apply(IndelsClassifier, axis=1)

    # DBS class
    indel_df1 = indel_df[indel_df['CLASS'] == 'DBS']
    indel_df1['CLASS'] = 'DBS'

    # INDEL class
    indel_df2 = indel_df[indel_df['CLASS'] != 'DBS']
    indel_df2['CLASS'] = 'INDEL'

    # merge both types
    indels = pd.concat([indel_df1, indel_df2, complex_indels_df, mnv_indels_df])

    indels.drop(['len_ref', 'len_alt'], axis=1, inplace=True)

    return indels


# main function
def format_hartwig(mutation_file, cnvs_file, purity_file, outfile):

    # load files and preformat them
    df, cnv_bed, purity_score, gender = load_files(mutation_file, cnvs_file, purity_file)

    # this is the sample column
    lastcol = list(df.columns)[-1]

    # get total reads
    df_reads = df.apply(get_reads, axis=1, args=([lastcol]))

    # select whether we have SNVs or others
    df_reads['len_alt'] = df_reads['ALT'].str.len()

    # number of characters in ref
    df_reads['len_ref'] = df_reads['REF'].str.len()

    # first classification between SNV and others
    df_reads['TYPE'] = df_reads.apply(
        lambda x: 'SNV' if ((x['len_alt'] == 1) and (x['len_ref'] == 1) and (x['ALT'] != '-') and (x['REF'] != '-')) else 'INDEL', axis=1
    )

    df_reads['pos-1'] = df_reads['POS'] - 1

    # get the triplet
    df_reads['TRIPLET'] = df_reads.apply(lambda x: hg19(x['CHROM'], x['pos-1'], 3), axis=1)
    df_reads['EXTENDED'] = df_reads.apply(lambda x: hg19(x['CHROM'], int(x['POS']) - 2, 5), axis=1)

    snv_df = df_reads[df_reads['TYPE'] != 'INDEL']
    snv_df['CLASS'] = 'SNV'
    snv_df['VARIANT_CLASS'] = snv_df.apply(create_snv_class, axis=1)

    # classify indels
    indel_df = df_reads[df_reads['TYPE'] == 'INDEL']
    indels = indels_classification(indel_df)
    columns = indels.columns

    df_reads_merged = pd.concat([snv_df, indels], sort=True)
    df_reads_merged = df_reads_merged[columns]

    # assing the name of the sample
    df_reads_merged['sample'] = lastcol

    # create bed file
    mut_bed = BedTool.from_dataframe(df_reads_merged[[
        'CHROM', 'pos-1', 'POS', 'ref_reads', 'var_reads', 'VAF', 'total_reads', 'REF',
        'ALT', 'sample', 'TYPE', 'CLASS', 'VARIANT_CLASS', 'TRIPLET', 'EXTENDED'
    ]])

    # Remove unmappable regions
    mapped = get_mappable_regions(mut_bed)

    # intersect with CN data
    out = mapped.intersect(cnv_bed, wao=True)

    # merge to dataframe
    merge = out.to_dataframe(names=[
        'CHROM', 'POS-1', 'POS', 'REF_COUNTS', 'VAR_COUNTS', 'VAF', 'TOTAL_READS', 'REF', 'ALT', 'SAMPLE', 'TYPE',
        'CLASS', 'VARIANT_CLASS', 'TRIPLET', 'EXTENDED', 'c1', 'p1', 'p2', 'MAJOR_CN_TEMP', 'actual_Baf', 'overlapp'
    ])

    # get the normal copy number values
    sex_chrom = ('Y', 'X')

    # get normal CN in the chromosome
    merge['NORMAL_CN'] = merge['CHROM'].apply(lambda x: 1 if x in sex_chrom and gender == "MALE" else 2)

    # add the purity score we got from PURPLE
    merge['PURITY'] = purity_score
    merge['GENDER'] = gender

    # get number of CNAs, if no overlapp then get the normal count
    merge['TOTAL_CN'] = merge.apply(get_major_cn, axis=1)

    # formula of allele specific copy number according to hartwig's people
    merge['MAJOR_CN'] = round(merge['actual_Baf'] * merge['TOTAL_CN']).astype(int)
    merge['MINOR_CN'] = round((1 - merge['actual_Baf']) * merge['TOTAL_CN']).astype(int)

    merge['CHROM'] = merge['CHROM'].apply(lambda x: 'chr{}'.format(x))

    # save files
    merge.dropna()[[
        'CHROM', 'POS', 'REF', 'ALT', 'TRIPLET', 'EXTENDED', 'CLASS', 'VARIANT_CLASS', 'SAMPLE', 'MAJOR_CN',
        'MINOR_CN', 'TOTAL_CN', 'NORMAL_CN', 'VAR_COUNTS', 'REF_COUNTS', 'GENDER', 'PURITY']].to_csv(
        outfile, sep='\t', index=False, header=True, compression='gzip'
    )

    # clean BedTools temp files
    pybedtools.cleanup()


if __name__ == '__main__':
    muts = sys.argv[1]
    cns = sys.argv[2]
    purity = sys.argv[3]
    outfile = sys.argv[4]

    format_hartwig(muts, cns, purity, outfile)
