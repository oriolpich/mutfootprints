# Import modules
import os
import sys
import gzip
from glob import glob

import seaborn as sns
import pandas as pd
import dill as pickle
import matplotlib as mpl
import matplotlib.pyplot as plt

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..',))
sys.path.insert(0, scripts_path)

from utils.plots_utils import config_params
from config import Conf


def list_MMR_genes():

    MMR_genes = 'data/repair_genes/repair_publication.tsv'
    MMR_df = pd.read_csv(MMR_genes, sep='\t', header=None)
    gene_list = MMR_df[0].tolist()

    return gene_list


def merge_coding_mutations_samples():

    toconcat = []
    path_annotated_samples = 'data/hartwig/samples/annotation/*.vcf.gz.annot.gz'

    for sample_file in glob(path_annotated_samples):
        sample_name = os.path.basename(sample_file).split('_')[1]
        df_sample = pd.read_csv(sample_file, sep='\t', low_memory=False)
        df_sample['SAMPLE'] = sample_name
        df_sample = df_sample[~df_sample['Func.wgEncodeGencodeBasicV19'].str.contains('ncRNA')]
        df_sample = df_sample[df_sample['Func.wgEncodeGencodeBasicV19'].str.contains('exonic|splicing')]
        toconcat.append(df_sample)

    coding_seq_samples = pd.concat(toconcat)

    return coding_seq_samples


def do_plot(
        exposures_pan, sig, MMR_notaffected_not_tzm_variants,
        MMR_affected_not_tzm_variants,
        MMR_notaffected_tzm_variants, MMR_affected_tzm_variants):

    config_params(6.5)

    fig, ax = plt.subplots(1, 1, figsize=(1.25, 1.5))
    sns.stripplot(
        data=[
            exposures_pan[MMR_notaffected_not_tzm_variants].loc[sig],
            exposures_pan[MMR_affected_not_tzm_variants].loc[sig],
            exposures_pan[MMR_notaffected_tzm_variants].loc[sig].tolist(),
            exposures_pan[MMR_affected_tzm_variants].loc[sig].tolist(), ],
        jitter=0.3, s=2, lw=0.5,
        color='#800080ff'
    )

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.xticks(
        [0, 1, 2, 3],
        [
            'TZM untreated, MMR-def (n={})'.format(len(MMR_notaffected_not_tzm_variants)),
            'TZM untreated, MMR-def (n={})'.format(len(MMR_affected_not_tzm_variants)),
            'TZM treated, no MMR-def (n={})'.format(len(MMR_notaffected_tzm_variants)),
            'TZM treated, MMR-def  (n={})'.format(len(MMR_affected_tzm_variants)),
        ],
        rotation=90
    )

    plt.ylabel('TMZ related SBS')
    plt.savefig('figures/TZM_treated.svg')
    plt.close()


def run():

    # signature similar to TMZ
    sig = '53_SBS11_0.988806_1'

    total_exonid = merge_coding_mutations_samples()
    MMR_list = list_MMR_genes()
    total_exonid['REPAIR'] = total_exonid['Gene.wgEncodeGencodeBasicV19'].apply(
        lambda x: 'MMR' if x in MMR_list else 'NO_MMR'
    )
    treated_samples = pickle.load(gzip.open(Conf['treatment_specific_drug']))
    exposures_path = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/'
    exposures_path += 'snvs/exposures/Pan_full/Pan_full.exposures.tsv'
    exposures_pan = pd.read_csv(exposures_path, sep='\t', index_col=0)

    # remove synonymous SNVs
    total_exonid = total_exonid[total_exonid['ExonicFunc.wgEncodeGencodeBasicV19'] != 'synonymous SNV']
    tzm_treated = [i for i in treated_samples['TEMOZOLOMIDE']['Pan']['YES'] if i in exposures_pan.columns]
    not_tzm_treated = [i for i in list(treated_samples['TEMOZOLOMIDE']['Pan']['NO']) if i in exposures_pan.columns]
    MMR_affected_not_tzm_variants = total_exonid[
        (total_exonid['SAMPLE'].isin(not_tzm_treated)) &
        (total_exonid['REPAIR'] == 'MMR')]['SAMPLE'].drop_duplicates()
    MMR_notaffected_not_tzm_variants = total_exonid[
        (total_exonid['SAMPLE'].isin(not_tzm_treated)) &
        (~total_exonid['SAMPLE'].isin(MMR_affected_not_tzm_variants))]['SAMPLE'].drop_duplicates()
    MMR_affected_tzm_variants = total_exonid[
        (total_exonid['SAMPLE'].isin(tzm_treated)) &
        (total_exonid['REPAIR'] == 'MMR')]['SAMPLE'].drop_duplicates()
    MMR_notaffected_tzm_variants = total_exonid[
        (total_exonid['SAMPLE'].isin(tzm_treated)) &
        (~total_exonid['SAMPLE'].isin(MMR_affected_tzm_variants))]['SAMPLE'].drop_duplicates()

    # generate the plot
    do_plot(
        exposures_pan,
        sig,
        MMR_notaffected_not_tzm_variants,
        MMR_affected_not_tzm_variants,
        MMR_notaffected_tzm_variants,
        MMR_affected_tzm_variants
    )


if __name__ == '__main__':
    run()
