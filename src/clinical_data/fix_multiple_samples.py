# Import modules
import os
import sys
from collections import defaultdict
from datetime import date

import pandas as pd

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf


# get the date
def create_date(x):
    if (str(x) != 'nan') and (str(x) != 'unknown'):
        date1 = date(int(x.split('-')[0]), int(x.split('-')[1]), int(x.split('-')[2]))
    else:
        date1 = 'unknown'
    return date1


# get path to treatments
def add_path_treatments(df, dic_updated_treatment):
    vals = dic_updated_treatment.get(df['sampleId'], df['preTreatments'])
    return vals


# time since we got the biopsy
def time_since_today(df):
    if df['BiopsyTime'] != 'unknown':
        distance = abs(df['BiopsyTime'] - df['todayDate']).days
    else:
        distance = 'unknown'

    return distance


def fix_treatments_multiple_samples():

    metadata_file = Conf['path_metadata']
    pat = pd.read_csv(metadata_file, sep='\t')

    # replace some bone soft tissues that where misspelled
    pat['primaryTumorLocation'] = pat['primaryTumorLocation'].replace('Bone/soft tissue', 'Bone/Soft tissue')
    # fix primary location to have compatible names
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation'].apply(
        lambda x: str(x).replace(' ', '-').replace('/', '-')
    )
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('Head-and-Neck', 'Head-and-neck')
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('nan', 'Unknown')
    pat['primaryTumorLocation_fixed'] = pat['primaryTumorLocation_fixed'].replace('CUP', 'Unknown')

    # create the dateTime from the biopsy date
    pat['BiopsyTime'] = pat['biopsyDate'].apply(create_date)

    # fix today's date
    pat['todayDate'] = date(2019, 5, 2)

    # date of treatment
    pat['treatmentStart'] = pat['treatmentStartDate'].apply(create_date)

    # date of end treatment
    pat['treatmentEnd'] = pat['treatmentEndDate'].apply(create_date)

    # days since the biopsy was taken
    pat['days_biopsy'] = pat.apply(time_since_today, axis=1)

    # We have to fix those patients without preTreatments but that they
    # have more than one sample and the first one is labelled as treated
    patient_occurrences = pat[~pat['biopsyDate'].isnull()]['patientId'].value_counts().to_dict()
    repeated_patient = []
    for p, times in patient_occurrences.items():
        if times > 1:
            repeated_patient.append(p)

    reppat = pat[pat['patientId'].isin(repeated_patient)]

    # remove those where biopsyDate is unknown since we cannot assess anything
    reppat = reppat[~reppat['biopsyDate'].isnull()]

    dic_updated_treatment = defaultdict()

    for patient, data in reppat.groupby(by='patientId'):

        # if we cannot guess which biopsy comes after the other, we cannot assess anything
        # But at least
        if 'unknown' not in data['days_biopsy'].tolist():

            # sort based on how old are the biopsies. This way the first biopsy is the first
            data.sort_values(by='days_biopsy', ascending=False, inplace=True)
            all_treat = []

            # if we have more than one biopsy per patient
            if len(data) > 1:
                for ix, (i, row) in enumerate(data.iterrows()):

                    # for the oldest one, we get this info
                    if ix == 0:

                        # this is what was given to them before the biopsy
                        treatment_old = row['preTreatments']

                        # this is what they plan to give them and should be annotated to subsequent samples
                        treatment_new = row['treatment']

                        # the preTreatment for the next biopsy should be both!
                        all_treat = [treatment_old, treatment_new]
                        all_treat = '/'.join([str(d) for d in all_treat if str(d) != 'nan'])

                        try:
                            # days since start of the treatment to today
                            start_treatment = abs(row['treatmentStart'] - row['todayDate']).days
                        except Exception:
                            start_treatment = 'unknown'

                    else:
                        # check if the time of treatment is greater than the time of the next biopsy,
                        # otherwise we cannot claim it might have had an impact on the biopsy
                        if start_treatment != 'unknown':

                            # if the biopsy was taken after the start of the treatment
                            if (start_treatment - row['days_biopsy']) > 0:
                                dic_updated_treatment[row['sampleId']] = all_treat

                                # update the treatments for the next iteration, if any
                                all_treat = [all_treat, row['treatment']]
                                all_treat = '/'.join([str(d) for d in all_treat if str(d) != 'nan'])

    pat['preTreatments_fixed'] = pat.apply(add_path_treatments, axis=1, args=(dic_updated_treatment,))
    pat.to_csv(Conf['preTreatments'], sep='\t', header=True, index=False)


if __name__ == '__main__':
    fix_treatments_multiple_samples()
