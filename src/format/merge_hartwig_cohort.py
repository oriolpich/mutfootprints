import pandas as pd
from glob import glob
from tqdm import tqdm
import sys
import numpy as np
from collections import defaultdict
import os
from scipy.stats import iqr

# insert path to root
scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))

sys.path.insert(0, scripts_path)
from utils.metadata import return_metadata


def merge_individual_samples_hartwig(path_hartwig, outpath):
    """
    This should merge all tumor types, remove those which belong to a primary,
    and then prepare for the extraction.
    :param path_hartwig:
    :return:
    """

    os.makedirs(outpath, exist_ok=True)

    # metadata of biopsy
    dic_primary_full, dic_secondary_fixed = return_metadata()

    all_hartwig = glob(path_hartwig)
    dic_muts = defaultdict(int)

    all_files = []
    for f in tqdm(all_hartwig, total=len(all_hartwig)):

        # remove those that do not belong to a good origin
        sample_name = os.path.basename(f).split('.')[0].split('_')[1]
        if sample_name in dic_primary_full:
            df = pd.read_csv(f, sep='\t')
            all_files.append(df)
            dic_muts[sample_name] = len(df)

    hart_df = pd.concat(all_files)

    print('sorting...')
    hart_df.sort_values(by=['CHROM', 'POS'], inplace=True)

    hart_df['PRIMARY'] = hart_df['SAMPLE'].map(dic_primary_full)
    hart_df['PRIMARY_EXTENDED'] = hart_df['SAMPLE'].map(dic_secondary_fixed)

    # we will split indels and SNVs, without considering hypermutants
    for ttype, data in tqdm(hart_df.groupby(by='PRIMARY')):

        outfile = '{}/{}.all.gz'.format(outpath, ttype)
        data.to_csv(outfile, sep='\t', header=True, index=False, compression='gzip')

        outfile = '{}/{}.snvs.gz'.format(outpath, ttype)
        data[data['CLASS'] == 'SNV'].to_csv(
            outfile, sep='\t', header=True, index=False, compression='gzip'
        )

        outfile = '{}/{}.indels.gz'.format(outpath, ttype)
        data[data['CLASS'] == 'INDEL'].to_csv(
            outfile, sep='\t', header=True, index=False, compression='gzip'
        )

        outfile = '{}/{}.dbs.gz'.format(outpath, ttype)
        data[data['CLASS'] == 'DBS'].to_csv(
            outfile, sep='\t', header=True, index=False, compression='gzip'
        )

        outfile = '{}/{}.mnvs.gz'.format(outpath, ttype)
        data[data['CLASS'] == 'MNV'].to_csv(
            outfile, sep='\t', header=True, index=False, compression='gzip'
        )

        outfile = '{}/{}.complex_indels.gz'.format(outpath, ttype)
        data[data['CLASS'] == 'COMPLEX_INDELS'].to_csv(
            outfile, sep='\t', header=True, index=False, compression='gzip'
        )

    # 2- save ALL file

    # Pan No Skin and no Unknown for SigProfiler extraction
    data = hart_df[
        (hart_df['PRIMARY'] != 'Skin') &
        (hart_df['PRIMARY'] != 'nan') &
        (hart_df['PRIMARY'] != 'Unknown') &
        (hart_df['PRIMARY'] != 'CUP')
        ]

    outfile = '{}/PanNoSkinNoUnknown.snvs.gz'.format(outpath)
    data[data['CLASS'] == 'SNV'].to_csv(
        outfile, sep='\t', header=True, index=False, compression='gzip'
    )

    outfile = '{}/PanNoSkinNoUnknown.indels.gz'.format(outpath)
    data[data['CLASS'] == 'INDEL'].to_csv(
        outfile, sep='\t', header=True, index=False, compression='gzip'
    )

    outfile = '{}/PanNoSkinNoUnknown.dbs.gz'.format(outpath)
    data[data['CLASS'] == 'DBS'].to_csv(
        outfile, sep='\t', header=True, index=False, compression='gzip'
    )

    outfile = '{}/PanNoSkinNoUnknown.mnvs.gz'.format(outpath)
    data[data['CLASS'] == 'MNV'].to_csv(
        outfile, sep='\t', header=True, index=False, compression='gzip'
    )

    outfile = '{}/PanNoSkinNoUnknown.complex_indels.gz'.format(outpath)
    data[data['CLASS'] == 'COMPLEX_INDELS'].to_csv(
        outfile, sep='\t', header=True, index=False, compression='gzip'
    )

    # define hypermutators based on SNV count. An hypermutated sample is a sample
    # with the total number of SNVs higher than 2.5IQR from the median of the entire distribution
    sample_counts = hart_df[hart_df['CLASS'] == 'SNV']['SAMPLE'].value_counts().to_dict()
    IQR = iqr(list(sample_counts.values()))
    median = np.median(list(sample_counts.values()))

    cutoff = median + 2.5 * IQR

    hypersamples = [s for s, v in sample_counts.items() if v > cutoff]
    nothypersamples = [s for s, v in sample_counts.items() if v <= cutoff]

    hypers = [hypersamples, nothypersamples]
    labels = ['Hypermutated', 'NotHypermutated']

    ttype = 'Pan'

    # Split into hypermutants and no hypermutants for Signature Analyzer extraction
    for h, l in zip(hypers, labels):

        subsdata = hart_df[hart_df['SAMPLE'].isin(h)]

        outfile = '{}/{}{}.all.gz'.format(outpath, ttype, l)
        subsdata.to_csv(
            outfile, sep='\t', header=True, index=False, compression='gzip'
        )

        outfile = '{}/{}{}.snvs.gz'.format(outpath, ttype, l)
        subsdata[subsdata['CLASS'] == 'SNV'].to_csv(
            outfile, sep='\t', header=True, index=False, compression='gzip'
        )

        outfile = '{}/{}{}.indels.gz'.format(outpath, ttype, l)
        subsdata[subsdata['CLASS'] == 'INDEL'].to_csv(
            outfile, sep='\t', header=True, index=False, compression='gzip'
        )

        outfile = '{}/{}{}.dbs.gz'.format(outpath, ttype, l)
        subsdata[subsdata['CLASS'] == 'DBS'].to_csv(
            outfile, sep='\t', header=True, index=False, compression='gzip'
        )

        outfile = '{}/{}{}.mnvs.gz'.format(outpath, ttype, l)
        subsdata[subsdata['CLASS'] == 'MNV'].to_csv(
            outfile, sep='\t', header=True, index=False, compression='gzip'
        )

        outfile = '{}/{}{}.complex_indels.gz'.format(outpath, ttype, l)
        subsdata[subsdata['CLASS'] == 'COMPLEX_INDELS'].to_csv(
            outfile, sep='\t', header=True, index=False, compression='gzip'
        )

    outfile = '{}/Pan.all.gz'.format(outpath)
    hart_df.to_csv(
        outfile, sep='\t', header=True, index=False, compression='gzip'
    )

    outfile = '{}/Pan.snvs.gz'.format(outpath)
    hart_df[hart_df['CLASS'] == 'SNV'].to_csv(
        outfile, sep='\t', header=True, index=False, compression='gzip'
    )

    outfile = '{}/Pan.indels.gz'.format(outpath)

    hart_df[hart_df['CLASS'] == 'INDEL'].to_csv(
        outfile, sep='\t', header=True, index=False, compression='gzip'
    )

    outfile = '{}/Pan.dbs.gz'.format(outpath)
    hart_df[hart_df['CLASS'] == 'DBS'].to_csv(
        outfile, sep='\t', header=True, index=False, compression='gzip'
    )

    outfile = '{}/Pan.mnvs.gz'.format(outpath)
    hart_df[hart_df['CLASS'] == 'MNV'].to_csv(
        outfile, sep='\t', header=True, index=False, compression='gzip'
    )

    outfile = '{}/Pan.complex_indels.gz'.format(outpath)
    hart_df[hart_df['CLASS'] == 'COMPLEX_INDELS'].to_csv(
        outfile, sep='\t', header=True, index=False, compression='gzip'
    )


if __name__ == '__main__':

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    merge_individual_samples_hartwig(input_path, output_path)
