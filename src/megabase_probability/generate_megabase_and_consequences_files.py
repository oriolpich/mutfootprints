import os
import sys

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf
import itertools
import subprocess
from pybedtools import BedTool
from bgvep.readers import BGPack
import pandas as pd
from collections import Counter, defaultdict
from bgreference import hg19
from tqdm import tqdm
import pickle, gzip
import numpy as np

def consequences_in_genes(genic_locations, df_mapp_bed, outfile_name):

    genic_full_regions = BedTool.from_dataframe(genic_locations[['chr', 'Gene start (bp)', 'Gene end (bp)']])

    # intersect mappable regions with all genic regions so that only regions with overlapping genes are considered,
    # to speed up calculations.
    all_genic_overlapp = df_mapp_bed.intersect(genic_full_regions, wa=True).to_dataframe(names=['chr', 'start', 'end', 'val',
                                                                                                'chr1', 'start1', 'end1', 'overlapp', 'ID', 'real_start'])

    conseq_type = {0: 'A',
                   1: 'C',
                   2: 'G',
                   3: 'T'}

    # set of consequence variants that we will consider as a protein-affecting
    consequence_wanted = {'start_lost', 'splice_region_variant', 'splice_donor_variant', 'stop_gained', 'stop_lost',
                              'missense_variant', 'splice_acceptor_variant'}

    # get the consequence type of our regions
    megabase_genic = defaultdict(dict)
    with BGPack('hg19', '88') as reader:
        for mb, data in tqdm(all_genic_overlapp.groupby(by='ID')):

            # counter of each type of mutation and how many times it affects the protein
            counter_feat = defaultdict(int)

            # go over each of the mappable intervals
            for i, row in data.iterrows():

                # for each of the positions within the intervals, check the consequence type
                for pos, cons in reader.get(row['chr'], row['real_start'], row['end']):

                    # each consequence type (most severe one), coming from A; C; G; T
                    for ix, c in enumerate(cons):

                        # if we are interested in the consequence
                        if c in consequence_wanted:
                            # keep the triplet
                            triplet = hg19(row['chr'], pos - 1, 3)
                            key = '{}_{}'.format(triplet, conseq_type[ix])
                            counter_feat[key] += 1

            # add the dict to each megabase
            megabase_genic[mb] = counter_feat

    pickle.dump(dict(megabase_genic), gzip.open('data/megabase_probability/{}.pckl.gz'.format(outfile_name),
                                                'wb'))

def calculate_possible_consequences():

    # mappable file
    mapp_file = 'data/megabase_probability/hg19.mappable.1Mb.windows.bed.extra.gz'

    # create bedfile
    df_mapp_bed = BedTool(mapp_file)

    # location of genes
    genic_locations = pd.read_csv('data/megabase_probability/gene_location_grch37.txt', sep='\t')
    genic_locations['chr'] = genic_locations['Chromosome/scaffold name'].apply(lambda x: 'chr{}'.format(x))

    # consequences in all genes
    consequences_in_genes(genic_locations, df_mapp_bed, "composition_affecting_genes_MB")

    # Cancer gene census
    CGC = Conf['CGC_path']
    dfcgc = pd.read_csv(CGC, sep ='\t')
    cgc_list = set(dfcgc['Gene Symbol'].tolist())
    genic_locations = genic_locations[genic_locations['Gene name'].isin(cgc_list)]

    # consequences in CGC
    consequences_in_genes(genic_locations, df_mapp_bed, "composition_affecting_CGCgenes_MB")


def get_mutations_MB_not_CGC():

    mapp_file = 'data/megabase_probability/hg19.mappable.1Mb.windows.bed.extra.gz'

    genes = pd.read_csv('data/megabase_probability/cgc_exonic_regions.tsv', sep ='\t', names = ['chr', 'p1', 'p2', 'strand', 'ID1', 'ID2', 'symbol'])

    genes['CHR'] = genes['chr'].apply(lambda x : 'chr{}'.format(x))

    genes_bed = BedTool.from_dataframe(genes[['CHR', 'p1', 'p2']])

    mapp_bedtool = BedTool(mapp_file)

    mapping_no_CGC = mapp_bedtool.subtract(genes_bed, sorted = True)

    mapp_no_CGC = mapping_no_CGC.to_dataframe(names= ['chr', 'start', 'end',  'val', 'chr1', 'start1','end1', 'overlapp', 'ID',
                                                   'real_start',])

    mapp_no_CGC.to_csv('data/megabase_probability/mappable_file.nocgc.bed.gz', sep ='\t', index = False,
                      header = False, compression = 'gzip')

def slicing_window(seq, n=3):

    it = iter(seq)
    result = ''.join(itertools.islice(it, n))

    if len(result) == n:
        yield result

    for elem in it:
        result = result[1:] + elem
        yield result

def counter_triplets_MB():

    # the file we just generated
    mapp_file = 'data/megabase_probability/hg19.mappable.1Mb.windows.bed.extra.gz'

    df_mapp = pd.read_csv(mapp_file, sep='\t', names=['chr', 'start', 'end', 'val', 'chr1', 'start1',
                                                      'end1', 'overlapp', 'ID', 'real_start'])

    df_mapp['len'] = df_mapp['end']-df_mapp['start']

    counter_per_megabase = defaultdict(dict)
    counter_nucl_per_megabase = defaultdict(int)

    for mb, region in tqdm(df_mapp.groupby(by='ID')):
        try:
            region['seq'] = region.apply(lambda x: hg19(x['chr'], x['start'],
                                                        x['end'] - x['start'] + 2), axis=1)
        except:
            region['seq'] = region.apply(lambda x: hg19(x['chr'], x['start'] + 1,
                                                        x['end'] - x['start']), axis=1)
        counter_region = Counter()
        for seq in region['seq'].tolist():
            sliced = Counter(list(slicing_window(seq)))
            counter_region += sliced

        counter_per_megabase[mb] = counter_region

        # count the length too
        counter_nucl_per_megabase[mb] = np.sum(region['len'].tolist())

    pickle.dump(dict(counter_per_megabase), gzip.open('data/megabase_probability/counter_1Mb.pckl.gz',
                                                      'wb'))
    pickle.dump(dict(counter_nucl_per_megabase), gzip.open('data/megabase_probability/mappable_counts_megabase_mutations.pckl.gz',
                                                      'wb'))
    total_count = defaultdict(int)
    for mb, d in counter_per_megabase.items():
        for triplet, c in d.items():
            total_count[triplet]+=c

    pickle.dump(dict(total_count), gzip.open('data/megabase_probability/counter_mappable.pckl.gz', 'wb'))

def run():

    print('counting stats per MB...')
    counter_triplets_MB()

    print('calculating possible consequences...')
    calculate_possible_consequences()

    print('removing mutations in CGS genes...')
    get_mutations_MB_not_CGC()

    cmd = '''zcat data/hartwig/origin/Pan.snvs.gz | tail -n +2 | awk '{OFS="\t"}{print $1, $2-1, $2, $3, $4, $8, $9, $18}' | intersectBed -a stdin -b data/megabase_probability/mappable_file.nocgc.bed.gz -wo | gzip > data/megabase_probability/pan.snvs.no_cgc_mappable_nottype.gz'''
    completed = subprocess.run(cmd, shell=True)

if __name__ == '__main__':
    run()
