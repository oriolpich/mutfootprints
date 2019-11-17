# Import modules
import os
import sys
from collections import defaultdict
import gzip
from itertools import chain
import pickle

import pandas as pd
import numpy as np

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf
from utils.serializable import leveled_default_dict
from drugs import composed_drugs, drugs_class_missing, fix_input_names


def load_FDA_classes():

    file_FDA = Conf['FDA']
    df = pd.read_csv(file_FDA, sep='\t', usecols=[
        'SUBSTANCENAME', 'NONPROPRIETARYNAME', 'PHARM_CLASSES'
    ])

    # remove duplicated lines
    df.drop_duplicates(inplace=True)

    df['CLASS'] = df['PHARM_CLASSES'].apply(lambda x: [t.split(' [EPC')[0] for t in str(x).split(',') if 'EPC' in t])
    df['CLASS'] = df['CLASS'].apply(lambda x: x if not isinstance(x, list) else x[0] if len(x) else np.nan)

    df = df[~df['CLASS'].isnull()]

    # FDA equivalence between drug name and drug class
    dic_t = dict(zip(df.SUBSTANCENAME, df.CLASS))

    # add the classes that were missing
    dic_t = {**drugs_class_missing, **dic_t}

    return dic_t


# get metadata
def load_metadata():

    file_metadata = Conf['preTreatments']

    pat = pd.read_csv(file_metadata, sep='\t')
    pat['primaryTumorLocation'] = pat['primaryTumorLocation'].replace('Bone/soft tissue', 'Bone/Soft tissue')

    # fix primary location
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation'].apply(
        lambda x: str(x).replace(' ', '-').replace('/', '-')
    )
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('Head-and-Neck', 'Head-and-neck')
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('nan', 'Unknown')
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('CUP', 'Unknown')

    dic_primary_full = dict(zip(pat.sampleId, pat.primaryTumorLocation_fixed))

    return pat, dic_primary_full


# fix correct names of drugs
def get_correct_names():

    metadata, dic_primary_full = load_metadata()
    dic_FDA = load_FDA_classes()

    # this will get the big class for each of the drugs
    drug_to_class = defaultdict(str)

    # this will store the drug given to the patient
    specific_drug = defaultdict(set)

    # fix each name to the canonical
    fixed_specific_name = defaultdict(set)

    # check not found items
    notfound = set()

    # dictionary where the key is the FDA class and the value the amount of patients treated per subtype
    typeDrug = leveled_default_dict(int, level=2)

    # analyze only those with pretreatments information
    treated = metadata[~metadata.preTreatments_fixed.isnull()]

    # each row is a sample
    for i, row in treated.iterrows():

        # split all treatments
        tspl = str(row['preTreatments_fixed']).split('/')
        sample = row['sampleId']

        # this list is erased for every sample
        treatment_done = []

        for t in tspl:

            # if it is in the composed dict, deconvolve it
            if t in composed_drugs:

                # loop over each of the drugs that belong to the formula
                for comp in composed_drugs[t]:

                    # make it upper and fix it if needed
                    tupper = fix_input_names.get(comp, comp).upper()

                    # check if if exists in FDA class
                    if tupper in dic_FDA:

                        specific_drug[sample].add(tupper)
                        drug_to_class[tupper] = dic_FDA[tupper]
                        fixed_specific_name[t].add(tupper)
                        class_drug = dic_FDA[tupper]
                        typeDrug[class_drug][tupper] += 1

                    else:
                        notfound.add(comp)

            # if it has no composite
            else:

                # make it upper and fix it if needed
                tupper = fix_input_names.get(t, t).upper()

                # check if if exists in FDA class
                if tupper in dic_FDA:

                    drug_to_class[tupper] = dic_FDA[tupper]
                    specific_drug[sample].add(tupper)
                    class_drug = dic_FDA[tupper]
                    typeDrug[class_drug][tupper] += 1
                    fixed_specific_name[t].add(tupper)

                    if tupper not in treatment_done:
                        typeDrug[class_drug][tupper] += 1
                        treatment_done.append(tupper)

                # try the hydrochloride formula in FDA
                else:
                    new_t = '{} HYDROCHLORIDE'.format(tupper)
                    if new_t in dic_FDA:
                        drug_to_class[tupper] = dic_FDA[new_t]
                        specific_drug[sample].add(new_t)
                        class_drug = dic_FDA[new_t]
                        typeDrug[class_drug][new_t] += 1
                        fixed_specific_name[t].add(new_t)

                    # try the acetate formula in FDA
                    else:
                        new_t = '{} ACETATE'.format(tupper)
                        if new_t in dic_FDA:
                            drug_to_class[tupper] = dic_FDA[new_t]
                            class_drug = dic_FDA[new_t]
                            specific_drug[sample].add(new_t)
                            typeDrug[class_drug][new_t] += 1
                            fixed_specific_name[t].add(new_t)

                        else:
                            notfound.add(t)

    # check if all drugs have an FDA assignment or not
    print('Drugs without assignment {}'.format(len(notfound)))
    if len(notfound) > 0:
        print(notfound)
        sys.exit()

    return drug_to_class, specific_drug, fixed_specific_name, typeDrug


# split samples into groups of treated/untreated
def samples_treatments_FDA_groups(drug_to_class, specific_drug, typeDrug):

    metadata, dic_primary_full = load_metadata()

    topo_II = [
        'DOXORUBICIN HYDROCHLORIDE', 'ETOPOSIDE', 'TENIPOSIDE',
        'MITOXANTRONE', 'ANTHRACYCLINES', 'EPIRUBICIN'
    ]

    # samples without systemic treatment
    samples_without_treatment = set(metadata[(metadata['hasSystemicPreTreatment'] == 'No') &
                                             (metadata['preTreatments_fixed'].isnull())]['sampleId'])
    # ---------------
    # FDA big group
    # ---------------

    # we will iterate over each FDA class drug
    samples_tract_FDA_dict = leveled_default_dict(set, level=3)

    for class_drug_to_test in set(drug_to_class.values()):

        # loop over all samples and their treatments
        for sample, sets in specific_drug.items():

            # flag to detect if the sample has any treatment that belongs to the class
            treat_c = False

            # add whether the treatment belongs to the class of FDA we are currently testing
            for drug in sets:
                if drug_to_class[drug] == class_drug_to_test:

                    # add sample to dictionary
                    samples_tract_FDA_dict[class_drug_to_test][dic_primary_full[sample]]['YES'].add(sample)
                    samples_tract_FDA_dict[class_drug_to_test]['Pan']['YES'].add(sample)

                    treat_c = True

            unknown_treat = False

            # add whether is not treated
            if treat_c is False:

                # check if any of the drugs given is unknown and remove the samples, just in case
                for drug in sets:
                    if 'UNKNO' in drug:
                        unknown_treat = True
                    elif 'UNKNO' in class_drug_to_test:
                        unknown_treat = True
                if unknown_treat is False:
                    samples_tract_FDA_dict[class_drug_to_test][dic_primary_full[sample]]['NO'].add(sample)
                    samples_tract_FDA_dict[class_drug_to_test]['Pan']['NO'].add(sample)

        for sample in samples_without_treatment:
            samples_tract_FDA_dict[class_drug_to_test][dic_primary_full[sample]]['NO'].add(sample)
            samples_tract_FDA_dict[class_drug_to_test]['Pan']['NO'].add(sample)

        # merge all classes of topoisomerases
        topo_class = 'TOPOII'
        for sample, sets in specific_drug.items():

            treat_c = False

            # add whether is treated
            for drug in sets:
                if drug in topo_II:
                    samples_tract_FDA_dict[topo_class][dic_primary_full[sample]]['YES'].add(sample)
                    samples_tract_FDA_dict[topo_class]['Pan']['YES'].add(sample)
                    treat_c = True

            unknown_treat = False

            # add whether is not treated
            if treat_c is False:
                for drug in sets:
                    if 'UNKNO' in drug:
                        unknown_treat = True
                    elif 'UNKNO' in class_drug_to_test:
                        unknown_treat = True
                if unknown_treat is False:
                    samples_tract_FDA_dict[topo_class][dic_primary_full[sample]]['NO'].add(sample)
                    samples_tract_FDA_dict[topo_class]['Pan']['NO'].add(sample)

    # ===
    # Specific drug
    # ===
    # we will iterate over each specific drug class
    samples_tract_specific_FDA = leveled_default_dict(set, level=3)
    for specific_drug_testing in list(set(chain.from_iterable(specific_drug.values()))):

        for sample, sets in specific_drug.items():

            # by default is a not treated sample
            treated_sample = False

            for drug in sets:

                # check if the drug is the same
                if drug == specific_drug_testing:
                    samples_tract_specific_FDA[specific_drug_testing][dic_primary_full[sample]]['YES'].add(sample)
                    samples_tract_specific_FDA[specific_drug_testing]['Pan']['YES'].add(sample)

                    treated_sample = True

            # if the sample has not been treated
            if treated_sample is False:

                unknown_treat = False

                for drug in sets:

                    if 'UNKNO' in drug:
                        unknown_treat = True
                    elif 'UNKNO' in drug_to_class[drug]:
                        unknown_treat = True

                if unknown_treat is False:
                    samples_tract_specific_FDA[specific_drug_testing][dic_primary_full[sample]]['NO'].add(sample)
                    samples_tract_specific_FDA[specific_drug_testing]['Pan']['NO'].add(sample)

        # complement this with the samples without annotation of treatment
        for sample in samples_without_treatment:
            samples_tract_specific_FDA[specific_drug_testing][dic_primary_full[sample]]['NO'].add(sample)
            samples_tract_specific_FDA[specific_drug_testing]['Pan']['NO'].add(sample)

    # add Radiation information
    for i, row in metadata.iterrows():

        if row['hasRadiotherapyPreTreatment'] == 'Yes':
            samples_tract_FDA_dict['RADIATION'][row['primaryTumorLocation_fixed']]['YES'].add(row['sampleId'])
            samples_tract_FDA_dict['RADIATION']['Pan']['YES'].add(row['sampleId'])
            samples_tract_specific_FDA['RADIATION'][row['primaryTumorLocation_fixed']]['YES'].add(row['sampleId'])
            samples_tract_specific_FDA['RADIATION']['Pan']['YES'].add(row['sampleId'])

        else:
            samples_tract_FDA_dict['RADIATION'][row['primaryTumorLocation_fixed']]['NO'].add(row['sampleId'])
            samples_tract_FDA_dict['RADIATION']['Pan']['NO'].add(row['sampleId'])
            samples_tract_specific_FDA['RADIATION'][row['primaryTumorLocation_fixed']]['NO'].add(row['sampleId'])
            samples_tract_specific_FDA['RADIATION']['Pan']['NO'].add(row['sampleId'])

    # Remove cases where the patient has been treated with a drug from a similar family
    stringent_filter = leveled_default_dict(set, level=3)

    # DBS inductors prior knowledge
    breaks = [
        'Unknown', 'Topoisomerase Inhibitor', 'Anthracycline Topoisomerase Inhibitor',
        'Alkylating Drug', 'Platinum-based Drug', 'Radioactive Diagnostic Agent', 'TOPOII'
    ]

    for drug in samples_tract_specific_FDA.keys():

        # select the FDA class for the drug
        class_drug = drug_to_class[drug]

        # create a list of samples treated with the same FDA type drug
        treated_with_same_type = []
        for drugs_of_type in typeDrug[class_drug].keys():
            if drugs_of_type != drug:
                treated_with_same_type += samples_tract_specific_FDA[drugs_of_type]['Pan']['YES']

        # DO SPECIAL CASE WITH RADIATION
        if drug == 'RADIATION':
            for class_drug in breaks:
                for drugs_of_type in typeDrug[class_drug].keys():
                    if drugs_of_type != drug:
                        treated_with_same_type += samples_tract_specific_FDA[drugs_of_type]['Pan']['YES']

        # select only those which are treated only with the specific type of drug as a positive and
        # those which are not treated with a similar drug as a negative
        not_treated = [s for s in samples_tract_specific_FDA[drug]['Pan']['NO'] if s not in treated_with_same_type]
        only_treated_specific_drug = [
            s for s in samples_tract_specific_FDA[drug]['Pan']['YES'] if s not in treated_with_same_type
        ]

        stringent_filter[drug]['Pan']['NO'] = set(not_treated)
        stringent_filter[drug]['Pan']['YES'] = set(only_treated_specific_drug)

    # merge 5-FU and capecitabine into 5-FU/CAPE
    for ttype in set(list(dic_primary_full.values())):

        samples_cape = samples_tract_specific_FDA['CAPECITABINE'][ttype]['YES']
        samples_fluoro = samples_tract_specific_FDA['FLUOROURACIL'][ttype]['YES']

        samples_NOT_cape = samples_tract_specific_FDA['CAPECITABINE'][ttype]['NO']
        samples_NOT_fluoro = samples_tract_specific_FDA['FLUOROURACIL'][ttype]['NO']

        total_treated = samples_fluoro | samples_cape
        total_NOT_treated = samples_NOT_cape & samples_NOT_fluoro

        samples_tract_specific_FDA['5-FU_CAPE'][ttype]['YES'] = total_treated
        samples_tract_specific_FDA['5-FU_CAPE'][ttype]['NO'] = total_NOT_treated

        for sample in total_treated:
            samples_tract_specific_FDA['5-FU_CAPE']['Pan']['YES'].add(sample)

        for sample in total_NOT_treated:
            samples_tract_specific_FDA['5-FU_CAPE']['Pan']['NO'].add(sample)

    return samples_tract_FDA_dict, samples_tract_specific_FDA, stringent_filter


# refactoring of some of the treatments, adding FDA categories and create list of treated/untreated per ttype
def refactoring_and_listing():

    drug_to_class, specific_drug, fixed_specific_name, typeDrug = get_correct_names()

    # dictionary from drug to FDA class
    pickle.dump(drug_to_class, gzip.open('data/clinical_data/hartwig_FDA_drug2class.pckl.gz', 'wb'))
    # dictionary from specific drug to canonical drug
    pickle.dump(fixed_specific_name, gzip.open('data/clinical_data/hartwig_fixed_name_treatments.pckl.gz', 'wb'))
    # dictionary with a counter for each type of drug
    pickle.dump(typeDrug, gzip.open('data/clinical_data/hartwig_typeDrug.pckl.gz', 'wb'))

    samples_tract_FDA_dict, samples_tract_specific_FDA, stringent_filter = samples_treatments_FDA_groups(
        drug_to_class, specific_drug, typeDrug
    )

    # dictionary with a list of treated/untreated for each FDA classes
    pickle.dump(samples_tract_FDA_dict, gzip.open(Conf['treatment_FDA_drug'], 'wb'))
    # dictionary with a list of treated/untreated for each of the specific drugs
    pickle.dump(samples_tract_specific_FDA, gzip.open(Conf['treatment_specific_drug'], 'wb'))
    # dictionary with a list of treated/untreated for each of the specific drugs with a stringent filer
    pickle.dump(stringent_filter, gzip.open(Conf['treatment_specific_drug_stringent'], 'wb'))


if __name__ == '__main__':
    refactoring_and_listing()
