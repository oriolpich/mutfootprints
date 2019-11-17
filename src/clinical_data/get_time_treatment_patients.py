# Import modules
import os
import sys
import copy
import gzip
import pickle
from datetime import date

import pandas as pd
import numpy as np

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf
from utils.serializable import leveled_default_dict


# get the date
def create_date(x):
    if (str(x) != 'nan') and (str(x) != 'unknown'):
        date1 = date(int(x.split('-')[0]), int(x.split('-')[1]), int(x.split('-')[2]))
    else:
        date1 = 'unknown'
    return date1

# get the length of the treatment (start-end)
def treatment_time(df):
    if df['endTreat'] != 'unknown' and df['startTreat'] != 'unknown':
        distance = abs(df['startTreat'] - df['endTreat']).days
    else:
        distance = 'unknown'
    return distance


# load refactored drug name
def load_name_drugs():
    fixed_specific_name = pickle.load(gzip.open('data/clinical_data/hartwig_fixed_name_treatments.pckl.gz'))
    return fixed_specific_name

# count length of exposures to treatments
def counter_exposures_treatment(pat, preDrugs):

    # get how many days in total the patient has been under the exposure of the drug
    counter_length = leveled_default_dict(list, level=2)

    # count how many times the patient has been treated with the specific drug
    counter_times = leveled_default_dict(int, level=2)

    # count the time since the start of the treatment
    counter_since_start_treatment = leveled_default_dict(list, level=2)

    # count the time since the treatment was last given
    counter_since_last_end_treatment = leveled_default_dict(list, level=2)

    fixed_specific_name = load_name_drugs()

    # loop over each of the samples
    for i, row in pat.iterrows():

        # patient Id
        patient = row['patientId']

        # sample Id
        sample = row['sampleId']

        # treatments received
        treats = str(row['preTreatments_fixed'])

        if treats != 'nan':

            # loop over all the treatments they got
            preTreat = row['preTreatments_fixed'].split('/')

            if len(preTreat):

                # loop over each of the treatments. We use set to avoid duplications
                for treat in set(preTreat):

                    # this means we have info on this drug
                    subs = preDrugs[(preDrugs['name'] == treat) & (preDrugs['#patientId'] == patient)]
                    if len(subs) > 0:

                        # we might have fixed some of the names (for example, combinations of drugs).
                        # This returns a list with the names
                        treat_fixed = fixed_specific_name[treat]
                        for t_fix in treat_fixed:

                            # for each of the times the patient received the drug
                            for j, row_med in subs.iterrows():

                                # this is the perfect case where we know biopsy timing, start and end of treatment
                                if (row_med['days_treatment'] != 'unknown') & (str(row['BiopsyTime']) != 'unknown'):

                                    # cases where the end of the treatment happens before the biopsy is taken
                                    if row_med['endTreat'] <= row['BiopsyTime']:

                                        # get the days of this specific treatment in this sample
                                        counter_length[sample][t_fix].append(row_med['days_treatment'])

                                        # times this person has been treated
                                        counter_times[sample][t_fix] += 1

                                        # time since it has been treated to the biopsy
                                        time_since_treatment = abs(row_med['startTreat'] - row['BiopsyTime']).days
                                        counter_since_start_treatment[sample][t_fix].append(time_since_treatment)

                                        # time since the treatment has finished (end treatment to biopsy)
                                        time_since_treatment_finish = abs(row_med['endTreat'] - row['BiopsyTime']).days
                                        counter_since_last_end_treatment[sample][t_fix].append(
                                            time_since_treatment_finish
                                        )
                                    else:
                                        # we can rescue those cases where the start
                                        # of the treatment is before the biopsy,
                                        if row_med['startTreat'] < row['BiopsyTime']:
                                            len_treat = abs(row_med['startTreat'] - row['BiopsyTime']).days
                                            time_since_treatment = len_treat

                                            counter_length[sample][t_fix].append(len_treat)
                                            counter_times[sample][t_fix] += 1
                                            counter_since_start_treatment[sample][t_fix].append(time_since_treatment)

                                            # this means there is no distance between the treatment and the biopsy,
                                            # the treatment is still being given to the patient
                                            counter_since_last_end_treatment[sample][t_fix].append(0)
                                        # this means the treatment has happened after the biopsy was taken
                                        else:
                                            continue
                                # if the info is not perfect, we cannot know the length of the treatment
                                else:
                                    if (str(row['BiopsyTime']) != 'unknown') & \
                                            (str(row_med['startTreat']) != 'unknown'):
                                        if row_med['startTreat'] < row['BiopsyTime']:
                                            time_since_treatment = abs(row_med['startTreat'] - row['BiopsyTime']).days
                                            counter_since_start_treatment[sample][t_fix].append(time_since_treatment)

                                    # in other cases we cannot do anything about it
                                    else:
                                        continue

                    for t_fix in treat_fixed:
                        if t_fix not in counter_length[sample]:
                            counter_length[sample][t_fix].append(np.nan)
                        if t_fix not in counter_times[sample]:
                            counter_times[sample][t_fix] += 1
                        if t_fix not in counter_since_start_treatment[sample]:
                            counter_since_start_treatment[sample][t_fix].append(np.nan)


    # add the merge of both drugs
    dic_treat_specific_copy = copy.deepcopy(counter_since_start_treatment)
    for sample, d_len in dic_treat_specific_copy.items():
        for drug, l in d_len.items():
            if drug in ['CAPECITABINE', 'FLUOROURACIL']:
                counter_since_start_treatment[sample]['5-FU_CAPE'].extend(l)

    counter_length_copy = copy.deepcopy(counter_length)
    for sample, d_len in counter_length_copy.items():
        for drug, l in d_len.items():
            if drug in ['CAPECITABINE', 'FLUOROURACIL']:
                counter_length[sample]['5-FU_CAPE'].extend(l)

    counter_times_copy = copy.deepcopy(counter_times)
    for sample, d_len in counter_times_copy.items():
        for drug, l in d_len.items():
            if drug in ['CAPECITABINE', 'FLUOROURACIL']:
                counter_times[sample]['5-FU_CAPE']+=l

    pickle.dump(counter_since_start_treatment, gzip.open(
        'data/clinical_data/dic_counter_since_start_treatment.pckl.gz', 'wb'
    ))
    pickle.dump(counter_length, gzip.open('data/clinical_data/dic_counter_total_days_treatment.pckl.gz', 'wb'))
    pickle.dump(counter_times, gzip.open('data/clinical_data/dic_counter_times_treatment.pckl.gz', 'wb'))
    preDrugs.to_csv('data/clinical_data/preDrugs_complete.tsv', sep='\t', index=False, header=True)


# get the timing of exposure to drugs
def timing_exposures_drugs():

    preDrugs_file = Conf['path_predrugs_original']
    preDrugs = pd.read_csv(preDrugs_file, sep='\t')

    # forbidden dates (their mistake, unsolvable)
    forb = ['1901-01-01', '1900-01-01']

    for fbd in forb:
        preDrugs.replace(fbd, np.nan, inplace=True)

    # replace one obvious mistake I found in one case
    preDrugs.replace('1015-04-03', '2015-04-03', inplace=True)

    # in the preDrugs file only those with annotated pretreatments are annotated, however
    # for the cases where we have two samples we can infer the treatment start for the second one

    preTreatments_file = Conf['preTreatments']
    preTreatments_df = pd.read_csv(preTreatments_file, sep='\t')

    preTreatments_df['BiopsyTime'] = preTreatments_df['biopsyDate'].apply(create_date)

    not_annotated = preTreatments_df[
        (preTreatments_df['preTreatments'].isnull()) &
        (preTreatments_df['treatmentStart'] != 'unknown')][
        ['patientId', 'treatmentStart', 'treatmentEnd', 'treatment']
    ]

    not_annotated['type'] = 'NaN'
    not_annotated['mechanism'] = 'NaN'
    not_annotated.columns = ['#patientId', 'startDate', 'endDate', 'name', 'type', 'mechanism']
    preDrugs_complete = pd.concat([preDrugs, not_annotated])
    preDrugs_complete['startTreat'] = preDrugs_complete['startDate'].apply(create_date)
    preDrugs_complete['endTreat'] = preDrugs_complete['endDate'].apply(create_date)
    preDrugs_complete['days_treatment'] = preDrugs_complete.apply(treatment_time, axis=1)

    counter_exposures_treatment(preTreatments_df, preDrugs_complete)


if __name__ == '__main__':
    timing_exposures_drugs()
