# Import modules
import os
import sys

import pandas as pd
from datetime import date
import matplotlib.pyplot as plt

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from utils.plots_utils import config_params


def time_since_today(df):
    if df['BiopsyTime'] != 'unknown':
        distance = abs(df['BiopsyTime'] - df['todayDate']).days
    else:
        distance = 'unknown'
    return distance


def time_since_today_drug(df):
    if df['treatmentStartDate_Time'] != 'unknown':
        distance = abs(df['treatmentStartDate_Time'] - df['todayDate']).days
    else:
        distance = 'unknown'
    return distance


def time_since_today_drug_end(df):
    if df['treatmentEndDate_Time'] != 'unknown':
        distance = abs(df['treatmentEndDate_Time'] - df['todayDate']).days
    else:
        distance = 'unknown'
    return distance


def time_since_today_treat(df):
    if df['startTreatTime'] != 'unknown':
        distance = abs(df['startTreatTime'] - df['todayDate']).days
    else:
        distance = 'unknown'
    return distance


def read_exposures(exposures_path):

    df = pd.read_csv(exposures_path, sep='\t', index_col=0)
    df = df.T
    samples_exp = df.index.tolist()
    return df, samples_exp


def create_date(x):
    if (str(x) != 'nan') & (str(x) != 'unknown'):
        date1 = date(int(x.split('-')[0]), int(x.split('-')[1]), int(x.split('-')[2]))
    else:
        date1 = 'unknown'
    return date1


def get_pre_post(drug, df_treat, dup):

    keep_patients = []
    stringent = []

    samples_not_treated = []
    samples_treated = []
    space_between_biopsies = []
    time_treated = []
    start_treatment = []

    for patient, data in dup.groupby(by='patientId'):
        data = data[data['days_biopsy'] != 'unknown']
        if len(data) > 1:
            data.sort_values(by='days_biopsy', ascending=False, inplace=True)
            first_sample = data.iloc[0]
            second_sample = data.iloc[1]

            if (drug in str(second_sample['preTreatments_fixed'])) & (first_sample['hasSystemicPreTreatment'] == 'No'):
                treat_date = df_treat[df_treat['#patientId'] == patient]
                treatment_between = treat_date[
                    (treat_date['days_since_treat'] <= first_sample['days_biopsy']) &
                    (treat_date['days_since_treat'] > second_sample['days_biopsy'])
                    ]

                found = False
                if len(treatment_between):
                    for i, row in treatment_between.iterrows():
                        if drug in row['name']:

                            if row['days_treatment'] != 'unknown':
                                keep_patients.append(patient)
                                samples_not_treated.append(first_sample['sampleId'])
                                samples_treated.append(second_sample['sampleId'])

                                # days until the biopsy
                                space_between_biopsies.append(
                                    first_sample['days_biopsy'] - second_sample['days_biopsy']
                                )

                                # time treated
                                time_treated.append(int(row['days_treatment']))

                                # time between the first biopsy and treatment
                                start_treatment.append(
                                    first_sample['days_biopsy'] - row['days_since_treat']
                                )

                                found = True

                    # if for whatever reasons the treatment is not found, use the big metadata
                    if found is False:

                        if first_sample['days_since_end_treatment'] != 'unknown':
                            keep_patients.append(patient)

                            samples_not_treated.append(first_sample['sampleId'])
                            samples_treated.append(second_sample['sampleId'])

                            # days until biopsy
                            space_between_biopsies.append(first_sample['days_biopsy'] - second_sample['days_biopsy'])

                            print(
                                first_sample['sampleId'],
                                second_sample['preTreatments_fixed'],
                                first_sample['preTreatments_fixed']
                            )  # time treated
                            # time treated
                            time_treated.append(
                                int(first_sample['days_since_start_treatment'] - first_sample['days_since_end_treatment'])
                            )

                            # time between the first biopsy and treatment
                            start_treatment.append(
                                first_sample['days_biopsy'] - first_sample['days_since_start_treatment']
                            )
                else:
                    if first_sample['days_since_end_treatment'] != 'unknown':

                        keep_patients.append(patient)
                        samples_not_treated.append(first_sample['sampleId'])
                        samples_treated.append(second_sample['sampleId'])

                        # days until biopsy
                        space_between_biopsies.append(first_sample['days_biopsy'] - second_sample['days_biopsy'])

                        print(
                            first_sample['sampleId'],
                            second_sample['preTreatments_fixed'],
                            first_sample['preTreatments_fixed']
                        )  # time treated
                        time_treated.append(
                            int(first_sample['days_since_start_treatment'] - first_sample['days_since_end_treatment'])
                        )

                        # time between the first biopsy and treatment
                        start_treatment.append(first_sample['days_biopsy'] - first_sample['days_since_start_treatment'])

    return (
        keep_patients, stringent, samples_not_treated, samples_treated,
        space_between_biopsies, time_treated, start_treatment
    )


def get_patients_two_points():
    preTreatments_info = 'data/clinical_data/preDrugs_complete.tsv'
    df_treat = pd.read_csv(preTreatments_info, sep='\t', low_memory=False)
    df_treat['todayDate'] = date(2019, 5, 2)
    df_treat['startTreatTime'] = df_treat['startTreat'].apply(create_date)
    df_treat['days_since_treat'] = df_treat.apply(time_since_today_treat, axis=1)

    pat = pd.read_csv('data/clinical_data/preTreatments_fixed.tsv', sep='\t')
    # create the dateTime from the biopsy date
    pat['BiopsyTime'] = pat['biopsyDate'].apply(create_date)
    pat['todayDate'] = date(2019, 5, 2)
    # days since the biopsy was taken
    pat['days_biopsy'] = pat.apply(time_since_today, axis=1)
    pat['treatmentStartDate_Time'] = pat['treatmentStartDate'].apply(create_date)
    pat['treatmentEndDate_Time'] = pat['treatmentEndDate'].apply(create_date)
    pat['days_since_start_treatment'] = pat.apply(time_since_today_drug, axis=1)
    pat['days_since_end_treatment'] = pat.apply(time_since_today_drug_end, axis=1)

    two_samples = {s for s, v in pat['patientId'].value_counts().to_dict().items() if v > 1}
    dup = pat[pat['patientId'].isin(two_samples)]

    return df_treat, dup


def analysis_two_timepoints(drug, signature_dic, exposures_path):

    config_params()

    lower_drug = drug.lower().capitalize()
    df_treat, dup = get_patients_two_points()
    (keep_patients, stringent, samples_not_treated, samples_treated,
     space_between_biopsies, time_treated, start_treatment) = get_pre_post(lower_drug, df_treat, dup)
    df_exp, samples_exposed = read_exposures(exposures_path)

    fig, ax = plt.subplots(1, 1, figsize=(2, 2))
    plt.ylabel('{}-related mutations'.format(lower_drug.capitalize()))
    plt.xlabel('Days')

    for ix in range(len(samples_not_treated)):
        if (samples_not_treated[ix] in samples_exposed) & (samples_treated[ix] in samples_exposed):

            exposure_not_treated = df_exp.loc[samples_not_treated[ix]][signature_dic[drug]].sum()
            exposure_treated = df_exp.loc[samples_treated[ix]][signature_dic[drug]].sum()

            plt.plot(
                [0, space_between_biopsies[ix]],
                [exposure_not_treated, exposure_treated],
                color='#d95f0e',
                ls='--'
            )

            slope = (exposure_treated - exposure_not_treated) / space_between_biopsies[ix]
            start_t = start_treatment[ix]
            end_t = start_treatment[ix] + time_treated[ix]

            plt.plot(
                [start_t, end_t],
                [exposure_not_treated + slope * start_t, exposure_not_treated + slope * end_t],
                color='darkred',
                lw=2
            )
            plt.scatter(space_between_biopsies[ix], exposure_treated)

    plt.savefig('figures/{}.pre_post.svg'.format(drug))
    plt.close()


def run():

    signature_dic = {
        'CAPECITABINE': ['31_SBS17b_0.968799_1'],
        'OXALIPLATIN': ['14_1', '37_1'],
    }

    exposures_path = 'data/hartwig/signatures/extraction/results/SignatureAnalyzer/'
    exposures_path += 'snvs/exposures/Pan_full/Pan_full.exposures.tsv'

    for drug in ['CAPECITABINE', 'OXALIPLATIN']:
        analysis_two_timepoints(drug, signature_dic, exposures_path)


if __name__ == '__main__':
    run()
