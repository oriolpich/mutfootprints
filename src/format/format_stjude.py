# Import modules
import os
import sys
from glob import glob

import pandas as pd
from pybedtools import BedTool

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf


def create_snv_class(df):

    pyr = ['C', 'T']
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    x = df['TRIPLET']

    if x[1] in pyr:
        out = '{}[{}>{}]{}'.format(x[0], x[1], df['ALT'], x[2])
    else:
        out = '{}[{}>{}]{}'.format(rev[x[2]], rev[x[1]], rev[df['ALT']], rev[x[0]])
    return out


def return_mappable_file():

    mappable_file = Conf['mappability_file']
    map_bed = BedTool(mappable_file)
    return map_bed


def merge_stjude_cohort():

    os.makedirs("data/STJUDE/format/", exist_ok=True)
    map_bed = return_mappable_file()
    treated_samples = [
        "SJOS001101_M1",
        "SJOS001101_M2",
        "SJOS001101_M3",
        "SJOS001101_M4",
        "SJOS001101_M5",
        "SJOS001101_M6",
        "SJOS001101_M7",
        "SJOS001101_M8",
        "SJOS001105_R1",
        "SJOS001105_D2",
        "SJOS001107_M2",
        "SJOS010_M",
        "SJOS010_D",
        "SJOS001_M"
    ]
    toconcat = []
    allf = glob(Conf['Path_STJUDE'] + '*.hg19.gz')

    for f in allf:
        sample_name = os.path.basename(f).split('.')[0]
        if sample_name in treated_samples:

            df = pd.read_csv(f, sep='\t')
            toconcat.append(df)

    df_concat = pd.concat(toconcat)
    df_concat['POSITION-1'] = df_concat['POSITION'] - 1
    wanted_order = ['CHROMOSOME', 'POSITION-1', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'DONOR_ID']
    df_concat = df_concat[wanted_order]

    df_bed = BedTool.from_dataframe(df_concat)
    intersected = df_bed.intersect(map_bed, wa=True)

    mappable_df = intersected.to_dataframe(names=['CHROM', 'POS-1', 'POS', 'REF', 'ALT', 'SAMPLE', 'DONOR_ID'])
    # select whether we have SNVs or others
    mappable_df['len_alt'] = mappable_df['ALT'].str.len()

    # number of characters in ref
    mappable_df['len_ref'] = mappable_df['REF'].str.len()
    mappable_df['TYPE'] = mappable_df.apply(lambda x: 'SNV' if (x['len_alt'] == 1) & (x['len_ref'] == 1) and (x['ALT']!='-') and (x['REF'] !='-') else 'INDEL', axis=1)
    mappable_df[mappable_df['TYPE'] == 'SNV'][['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT']].to_csv(
        'data/STJUDE/format/STJUDE.to_deconstruct.tsv', sep='\t', index=False, header=False
    )

if __name__ == '__main__':
    merge_stjude_cohort()
