# Import modules
import os
import sys
from glob import glob

import click
import pandas as pd
from tqdm import tqdm

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from utils.metadata import return_metadata


# Global variables
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input-file', required=True, help='Input file', type=click.Path())
@click.option('-o', '--outpath', required=True, help='Output directory with the results', type=click.Path())
def merge_samples_timing(input_file, outpath):

    os.makedirs(outpath, exist_ok=True)
    dic_primary_full, _ = return_metadata()
    all_hartwig = glob(input_file)
    to_remove = [
        'EXTENDED', 'MAJOR_CN', 'MINOR_CN', 'TOTAL_CN', 'NORMAL_CN',
        'VAR_COUNTS', 'REF_COUNTS', 'GENDER', 'PURITY'
    ]

    all_files = []
    for f in tqdm(all_hartwig, total=len(all_hartwig)):

        # remove those that do not belong to a good origin
        sample_name = os.path.basename(f).split('_')[1].split('.')[0]
        if sample_name in dic_primary_full:
            df = pd.read_csv(f, sep='\t', low_memory=False)
            df.drop(to_remove, axis=1, inplace=True)
            all_files.append(df)

    hart_df = pd.concat(all_files)

    hart_df['PRIMARY'] = hart_df['SAMPLE'].map(dic_primary_full)
    hart_df.sort_values(by=['CHROM', 'POS'], inplace=True)

    outfile = '{}/Pan.snvs.gz'.format(outpath)
    hart_df[hart_df['CLASS'] == 'SNV'].to_csv(outfile, sep='\t', header=True, index=False, compression='gzip')

    outfile = '{}/Pan.indels.gz'.format(outpath)

    hart_df[hart_df['CLASS'] == 'INDEL'].to_csv(outfile, sep='\t', header=True, index=False, compression='gzip')

    outfile = '{}/Pan.dbs.gz'.format(outpath)
    hart_df[hart_df['CLASS'] == 'DBS'].to_csv(outfile, sep='\t', header=True, index=False, compression='gzip')

    # save ALL file
    outfile = '{}/Pan.all.gz'.format(outpath)
    hart_df.to_csv(outfile, sep='\t', header=True, index=False, compression='gzip')


if __name__ == '__main__':
    merge_samples_timing()
