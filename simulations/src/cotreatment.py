"""
Explain each signature as a mixture of concomitant treatments.
For each signature and a given collection of concomitant treatments,
we deem "efficiency" the capacity for a treatment to generate exposure of a specific signature.
"""

import os
import tqdm
import click
import numpy as np
import pandas as pd
import scipy.optimize
import dill as pickle
import gzip

from non_negative_least_squares import nnls_input
from recovery import signature_recovery


class EfficiencyDrugSpecific:

    def __init__(self, exposures_path, matrix_treatments_path, signature_index):

        self.exposures = pd.read_csv(exposures_path, sep='\t')
        self.exposures = self.exposures.transpose()
        self.matrix_treatments = pd.read_csv(matrix_treatments_path, sep='\t', index_col=0)
        self.signature_index = signature_index

    def parse_exposure(self):

        return self.exposures.iloc[:, self.signature_index]

    def treatfit(self, sig, placebo):

        # exposures dataframe
        exposure = self.parse_exposure()
        treatments = self.matrix_treatments[[sig, placebo]]

        # NNLS: transform the input
        E = exposure.values
        T = treatments.values

        A, b = nnls_input(T, E)

        # NNLS: custom bounds -- some values must be zero
        lb = np.zeros(A.shape[1])  # lower bound
        ub = lb + np.inf

        # NNLS: solver
        res = scipy.optimize.lsq_linear(A, b, bounds=(lb, ub))

        # collapse near-zero coefficients
        return dict(zip([sig, placebo], list(res.x)))

    def treatfit_many_at_once(self, sig, competitors):

        # exposures dataframe
        exposure = self.parse_exposure()
        treatments = self.matrix_treatments[[sig] + competitors]

        # NNLS: transform the input
        E = exposure.values
        T = treatments.values

        A, b = nnls_input(T, E)

        # NNLS: custom bounds -- some values must be zero
        lb = np.zeros(A.shape[1])  # lower bound
        ub = lb + np.inf

        # NNLS: solver
        res = scipy.optimize.lsq_linear(A, b, bounds=(lb, ub))

        # collapse near-zero coefficients
        return dict(zip([sig] + competitors, list(res.x)))


def origin(treatment, signature, exposures_path, discriminate_path, matrix_treatment_path):

    # show progress during execution
    print(treatment, signature)

    # determine efficiencies of competing matrix_treatments
    e = EfficiencyDrugSpecific(exposures_path, discriminate_path, matrix_treatment_path)

    # corner case 1
    if signature not in e.exposures.columns:
        return None, None

    # corner case 2
    if treatment not in e.matrix_treatments.columns:
        return None, None

    dftreatpool = e.select_pools(treatment, signature)

    # corner case 3
    if len(dftreatpool) == 0:
        return None, None

    x = e.treatfit(dftreatpool, [signature])
    x = dict(zip(list(x.index), list(x.iloc[:, 0].values)))

    # most likely origin
    origin = sorted(x.keys(), key=lambda k: x[k], reverse=True)[0]

    # efficiency rate between origin and reference treatment
    if (origin != treatment) and (x[treatment] == 0):
        score = np.inf
    else:
        score = x[origin] / x[treatment]
    return origin, score


def dump_efficiency_dict(exposures_path, discriminate_path, matrix_treatment_path, output):

    eff = EfficiencyDrugSpecific(exposures_path, discriminate_path, matrix_treatment_path)
    with gzip.open(output, 'wb') as fn:
        pickle.dump(eff, fn)


def create_random_matrix_treatments(sig, n_treated, replicate, overlap, exposure_folder, output_folder, symmetric=True):

    exposure_path = os.path.join(exposure_folder,
                                 f'{sig}.{n_treated}.{replicate}',
                                 f'{sig}.{n_treated}.{replicate}.exposures.tsv')
    df = pd.read_csv(exposure_path, sep='\t')
    dg = df.transpose()
    dg[sig] = dg.apply(lambda x: 1 if x.name.endswith('_t') else 0, axis=1)
    index_treated = dg.loc[dg.index.str.endswith('_t')].index.to_list()
    index_nontreated = dg.loc[~dg.index.str.endswith('_t')].index.to_list()
    for i in range(10):
        placebo_treated = np.random.choice(index_treated, size=len(index_treated) * overlap // 100, replace=False)
        if symmetric:
            placebo_nontreated = np.random.choice(index_nontreated, size=len(index_nontreated) * overlap // 100, replace=False)
        else:
            placebo_nontreated = np.random.choice(index_nontreated, size=n_treated - len(placebo_treated), replace=False)
        placebo = list(placebo_treated) + list(placebo_nontreated)
        dg[f'placebo_{i+1}'] = dg.apply(lambda x: int(x.name in placebo), axis=1)
    dg = dg[[sig] + [f'placebo_{i+1}' for i in range(10)]]
    symmetric_label = 'symmetric' if symmetric else 'asymmetric'
    output_path = os.path.join(output_folder,
                               f'{sig}.{n_treated}.{replicate}.{overlap}.{symmetric_label}.treatment_matrix.tsv')
    dg.to_csv(output_path, sep='\t', index=True)


def attribution(sig, n_treated, replicate, deconstruction_folder, matrix_treatment_folder, overlap, signatures, symmetric):

    symmetric_label = 'symmetric' if symmetric else 'asymmetric'
    exposures_path = os.path.join(deconstruction_folder,
                                  f'{sig}.{n_treated}.{replicate}',
                                  f'{sig}.{n_treated}.{replicate}.exposures.tsv')
    matrix_treatment_path = os.path.join(matrix_treatment_folder,
                                         f'{sig}.{n_treated}.{replicate}.{overlap}.{symmetric_label}.treatment_matrix.tsv')

    sig_index, cosine, twin_signals = signature_recovery(sig, n_treated, replicate, signatures, deconstruction_folder)

    # determine efficiencies of competing matrix_treatments

    e = EfficiencyDrugSpecific(exposures_path, matrix_treatment_path, sig_index)
    res_dict = {sig: [], 'placebo': []}

    for i in range(10):
        fit = e.treatfit(sig, f'placebo_{i+1}')
        res_dict[sig].append(fit[sig])
        res_dict['placebo'].append(fit[f'placebo_{i+1}'])
    return res_dict


def attribution_many_to_one(sig, n_treated, replicate, deconstruction_folder, matrix_treatment_folder, overlap, signatures, symmetric):

    symmetric_label = 'symmetric' if symmetric else 'asymmetric'
    exposures_path = os.path.join(deconstruction_folder,
                                  f'{sig}.{n_treated}.{replicate}',
                                  f'{sig}.{n_treated}.{replicate}.exposures.tsv')
    matrix_treatment_path = os.path.join(matrix_treatment_folder,
                                         f'{sig}.{n_treated}.{replicate}.{overlap}.{symmetric_label}.treatment_matrix.tsv')
    sig_index, cosine, twin_signals = signature_recovery(sig, n_treated, replicate, signatures, deconstruction_folder)

    # determine efficiencies of competing matrix_treatments

    e = EfficiencyDrugSpecific(exposures_path, matrix_treatment_path, sig_index)
    competitors = pd.read_csv(matrix_treatment_path, sep='\t', index_col=None)
    competitors = list(competitors.columns[2:])
    res_dict = {sig: [], 'placebo': []}

    for i, c in enumerate(competitors):
        fit = e.treatfit_many_at_once(sig, competitors[: i+1])
        res_dict[sig].append(fit[sig])
        del fit[sig]
        max = np.max(list(fit.values()))
        res_dict['placebo'].append(max)

    return res_dict


@click.group()
def cli():
    pass



@cli.command()
@click.option('--deconstruction_folder', type=click.Path())
@click.option('--treatment_matrix_folder', type=click.Path())
@click.option('--symmetric', is_flag=True)
def run_create_random_matrix_treatments(deconstruction_folder, treatment_matrix_folder, symmetric):

    for sig in tqdm.tqdm(['SBS9', 'SBS31']):
        for n_treated in tqdm.tqdm([10, 25, 50, 100, 150]):
            for replicate in range(1, 11):
                for overlap in [25, 50, 75]:
                    try:
                        create_random_matrix_treatments(sig, n_treated, replicate, overlap,
                                                        deconstruction_folder, treatment_matrix_folder,
                                                        symmetric=symmetric)
                    except FileNotFoundError:
                        print(sig, n_treated, replicate)


if __name__ == '__main__':

    cli()
