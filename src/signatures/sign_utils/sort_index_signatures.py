# Import modules
import os
import sys

import pandas as pd

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf


# read  canonical SBS signatures from PCAWG
def pcawg_canonical_snvs():

    pcawg_sbs_file = Conf['SBS_PCAWG']
    pcawg_snvs = pd.read_csv(pcawg_sbs_file)
    pcawg_snvs.index = pcawg_snvs.apply(
        lambda x: '{}[{}>{}]{}'.format(
            x['SubType'][0], x['SubType'][1], x['Type'][-1], x['SubType'][2]
        ), axis=1
    )

    pcawg_snvs.sort_index(inplace=True)
    pcawg_snvs.drop(['SubType', 'Type'], axis=1, inplace=True)
    pcawg_snvs = pcawg_snvs/pcawg_snvs.sum()
    pcawg_snvs.to_csv(pcawg_sbs_file.replace('.csv', '.indexed.csv'), sep='\t', index=False, header=True)

    pcawg_snvs = pd.read_csv(pcawg_sbs_file)
    # for comparisons with signatureanalyzer
    pcawg_snvs.index = pcawg_snvs.apply(
        lambda x: '{}{}{}{}'.format(
            x['SubType'][0], x['SubType'][1], x['SubType'][2], x['Type'][-1]
        ), axis=1
    )

    pcawg_snvs.sort_index(inplace=True)
    pcawg_snvs.drop(['SubType', 'Type'], axis=1, inplace=True)
    pcawg_snvs = pcawg_snvs/pcawg_snvs.sum()
    pcawg_snvs.to_csv(pcawg_sbs_file.replace('.csv', '.extended_indexed.csv'), sep='\t', index=False, header=True)


# the code might be redundant
def create_index_1536():

    pcawg_sbs_file = Conf['SBS_PCAWG']
    outfile = '{}/lego1536.index.tsv'.format(os.path.dirname(pcawg_sbs_file))
    nuc = ['A', 'C', 'G', 'T']
    midd = ['C', 'T']
    total = []

    # this is the reference
    for ref in midd:
        for mut in nuc:
            if mut != ref:
                for f in nuc:
                    for s in nuc:
                        for four in nuc:
                            for last in nuc:
                                a = '{}{}{}{}{}{}'.format(f, s, ref, four, last, mut)
                                total.append(a)

    with open(outfile, 'wt') as outf:
        for element in sorted(total):
            out = '{}\n'.format(element)
            outf.write(out)

    total = []
    # this is the reference
    for ref in midd:
        for mut in nuc:
            if mut != ref:
                for s in nuc:
                    for last in nuc:
                        a = '{}{}{}{}'.format(s, ref, last, mut)
                        total.append(a)

    outfile = '{}/lego96.index.tsv'.format(os.path.dirname(pcawg_sbs_file))
    with open(outfile, 'wt') as outf:
        for ix, element in enumerate(sorted(total)):
            out = '{}\t{}\n'.format(ix + 1, element)
            outf.write(out)


# for the GBM analysis
def create_signatures_for_deconstructsigs():
    pcawg_sbs_file = Conf['SBS_PCAWG']
    outfile = '{}/SigProfiler_SBS_signatures_2018_03_28.deconstructsigs.tsv'.format(os.path.dirname(pcawg_sbs_file))
    pcawg_snvs = pd.read_csv(pcawg_sbs_file)
    pcawg_snvs.index = pcawg_snvs.apply(
        lambda x: '{}[{}>{}]{}'.format(
            x['SubType'][0], x['SubType'][1], x['Type'][-1], x['SubType'][2]
        ), axis=1
    )
    pcawg_snvs.drop(['Type', 'SubType'], axis=1).T.to_csv(outfile, index=True, header=True, sep='\t')


def run():
    create_index_1536()
    pcawg_canonical_snvs()
    create_signatures_for_deconstructsigs()


if __name__ == '__main__':
    run()
