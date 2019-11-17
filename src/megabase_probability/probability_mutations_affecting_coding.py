import pandas as pd
import gzip
from collections import defaultdict
import numpy as np
from tqdm import tqdm
import pickle

def return_signature_probability(path_to_process, counter_mappable, signature):

    exp = pd.read_csv( path_to_process, sep='\t', index_col=0)

    vals = defaultdict(float)
    reverse = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    for i, row in exp.iterrows():

        key = '{}{}{}_{}'.format(i[0], i[2], i[-1], i[4])
        triplet = '{}{}{}'.format(i[0], i[2], i[-1])

        # reverse complementary
        rev_triplet = ''.join([reverse[ix] for ix in triplet[::-1]])
        key_rev = '{}_{}'.format(rev_triplet, reverse[i[4]])

        vals[key] += row[signature]/counter_mappable[triplet]
        vals[key_rev] += row[signature]/counter_mappable[rev_triplet]

    vals_norm = defaultdict(float)
    total_c = np.sum(list(vals.values()))
    for trip, c in vals.items():
        vals_norm[trip] = c / total_c

    return vals_norm


def return_mut_risk(path_to_process, treatment, ttype, signature, number_of_mutations):

    # we should normalize the MB
    mut_rate_normalized = defaultdict(float)
    real_mb_prob = mut_rate_ttype[treatment][ttype]

    total_in_mb = np.sum(list(mut_rate_ttype[treatment][ttype].values()))

    for mb, val in real_mb_prob.items():
        mut_rate_normalized[mb] = val / total_in_mb

    final_count_genic = 0
    final_count_cgc = 0

    signature_probability = return_signature_probability(path_to_process, counter_mappable, signature)

    # this is for genic
    for mb in tqdm(mut_rate_normalized.keys()):

        if mb in megabase_genic:

            # given the counts per megabase and how likely is the context of the signature,
            # calculate how likely is to observe a mutation in a context given the composition of the MB
            probability_context_mb = defaultdict(float)
            count_context_triplet = defaultdict(float)
            for context, value in counter_per_megabase[mb].items():
                for nucleotide in 'ACGT':
                    if nucleotide != context[1]:
                        mutation = '{}_{}'.format(context, nucleotide)
                        probability_context_mb[mutation] = value * signature_probability[mutation]

                        # keep the "fake" count of the mutation to occur
                        count_context_triplet[mutation] = value

            # normalize the probability. This means, if we were to have one mutation,
            # in which context they would be more likely to fall.
            probability_context_mb_norm = defaultdict(float)
            total_contex_mb = np.sum(list(probability_context_mb.values()))
            for trip, c in probability_context_mb.items():
                probability_context_mb_norm[trip] = c / total_contex_mb

            # we calculate the proportion of sites of the context that are affecting (given our signature), that is, the
            # count of the genic by the total in the megabase
            probabilities_in_mb = defaultdict(float)
            for trip, c in megabase_genic[mb].items():
                probabilities_in_mb[trip] = c / count_context_triplet[trip]

            # finally, to know exactly the probability, assuming that one mutation is falling in the MB, we should multiply
            # the chances of the mutation falling in a particular context by the probability of this context affecting a gene
            probability_mut_coding = 0
            for triplet, likely_coding in probabilities_in_mb.items():
                probability_mut_coding += likely_coding * probability_context_mb_norm[triplet]

            # of course, this should be normalized by the probability of the MB to be mutated and by the number of mutations observed
            norm_prob_in_MB = probability_mut_coding * mut_rate_normalized[mb] * number_of_mutations
            final_count_genic += norm_prob_in_MB

            # do the same for cancer gene census
            probabilities_in_mb_cgc = defaultdict(float)
            if mb in megabase_CGC_genic:

                for trip, c in megabase_CGC_genic[mb].items():
                    probabilities_in_mb_cgc[trip] = c / count_context_triplet[trip]

                probability_mut_coding = 0
                for triplet, likely_coding in probabilities_in_mb_cgc.items():
                    probability_mut_coding += likely_coding * probability_context_mb_norm[triplet]

                norm_prob_in_MB_CGC = probability_mut_coding * mut_rate_normalized[mb] * number_of_mutations
                final_count_cgc += norm_prob_in_MB_CGC

    return final_count_genic, final_count_cgc


# files needed

# counter of triplets per megabase
counter_per_megabase = pickle.load(gzip.open('data/megabase_probability/counter_1Mb.pckl.gz'))

# counter of triplets mappable
counter_mappable = pickle.load(gzip.open('data/megabase_probability/counter_mappable.pckl.gz'))

# number of triplets affecting genes per megabase
megabase_genic = pickle.load(gzip.open('data/megabase_probability/composition_affecting_genes_MB.pckl.gz', ))

# normalized mutation rate per megabase (per tumor type)
mut_rate_ttype = pickle.load(gzip.open('data/megabase_probability/normalized_megabase_mutations_ttype.pickle.gz'))

# number of triplets affecting genes per megabase
megabase_CGC_genic = pickle.load(gzip.open('data/megabase_probability/composition_affecting_CGCgenes_MB.pckl.gz', ))

# path where processes files are located
path_to_process = "data/hartwig/signatures/extraction/results/SigProfiler/snvs/processes/PanNoSkinNoUnknown.snvs/PanNoSkinNoUnknown.snvs.processes.tsv"
wanted_ttypes = ['Ovary', 'Colon-Rectum', 'Esophagus', 'Lung', 'Urinary-tract', 'Uterus',   'Breast', 'NET']

sigprofiler_sigs = {'CAPECITABINE': '19_SBS17b_0.961548_0.99',
                    'CARBOPLATIN': '1_SBS31_0.968153_0.98',
                    'CISPLATIN': '1_SBS31_0.968153_0.98',
                    'OXALIPLATIN': '20_0.92',
                    'AGING': '23_SBS1_0.996494_1.0',
                    'SMOKING': '17_SBS4_0.935484_0.97'}

number_of_mutations = 1
res_ttype = defaultdict(lambda: defaultdict(dict))

for ttype in wanted_ttypes:

    for treatment in tqdm(sigprofiler_sigs.keys()):
        gene, cgc = return_mut_risk(path_to_process, treatment, ttype, sigprofiler_sigs[treatment], 1)
        res_ttype[ttype][treatment]['gene'] = gene
        res_ttype[ttype][treatment]['cgc'] = cgc

pickle.dump(dict(res_ttype), gzip.open('data/megabase_probability/results_probability_sigpro_20190502.pckl.gz', 'wb'))
