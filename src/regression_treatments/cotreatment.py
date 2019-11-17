"""
Explains each signature as a mixture of co-administered treatments.
For each signature and the collection of co-administered treatments that are related to it,
we deem "efficiency" the capacity for a treatment to generate exposure of a specific signature.
"""

import os
import sys

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]


scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

import click
import numpy as np
import pandas as pd
import scipy.optimize
import dill as pickle
import gzip

from non_negative_least_squares import nnls_input, nnls_input_bootstrap
from utils.metadata import return_metadata


# global namespace
PVALUE = 0.001
EFFECT_SIZE = 2
OVERLAP = 0.4


def lowerize(string):

    return '_'.join(list(map(str.lower, string.split(' '))))


def upperize(string):

    return ' '.join(list(map(str.upper, string.split('_'))))


class EfficiencyDrugSpecific:

    def __init__(self, exposures_path, discriminate_path, matrix_treatments_path, ttype='Pan'):

        self.exposures = pd.read_csv(exposures_path, sep='\t', index_col=0)
        self.exposures = self.exposures.transpose()
        self.discriminate = pd.read_csv(discriminate_path, sep='\t')
        self.discriminate = self.discriminate[self.discriminate['ttype'] == ttype]
        self.discriminate = self.discriminate.astype({'signature': str, 'effect_size': float, 'pvals': float,
                                                      'treatment': str, 'ttype': str, 'total_treated': float})
        self.treatments = pd.read_csv(matrix_treatments_path, sep='\t', index_col=0)

        self.ttype = ttype

        dic_ttypes, _ = return_metadata()

        if ttype == 'Pan':
            self.samples = dic_ttypes.keys()
        else:
            self.samples = [k for k, v in dic_ttypes.items() if v == ttype]

        self.effect_thresh = EFFECT_SIZE   # min effect size
        self.pvalue_thresh = PVALUE        # empirical_pvalue threshold
        self.overlap_thresh = OVERLAP      # min overlap to consider cotreatments

    def run(self):
        # TODO: this is not used

        df = self.discriminate
        df = df[(df.effect_size > self.effect_thresh) & (df.pvals < self.pvalue_thresh) & (df.ttype == self.ttype)]

        efficiency_dict = {}
        for t in df.treatment.unique():
            for s in df[df.treatment == t]['signature'].unique():
                treatpool = self.select_pools(t, s)
                x = self.treatfit(treatpool, [s])
                efficiency_dict[(t, s)] = dict(zip(list(map(lowerize, x.index)), x.iloc[:, 0].values))
        return efficiency_dict

    def neighbours(self, treatment):

        neighbours_dict = {}
        df = self.treatments
        df = df[df[treatment] == 1]
        for col in df.columns:
            if len(df) == 0:
                overlap = 0
            else:
                overlap = df[col].values.sum() / len(df)
            if overlap > 0:
                neighbours_dict[col] = overlap
        return [k for k, v in neighbours_dict.items() if v > self.overlap_thresh]

    def select_pools(self, treatment, signature):

        treatpool = self.neighbours(treatment)
        df = self.discriminate.copy()
        df = df[(df.effect_size > self.effect_thresh) & (df.pvals < self.pvalue_thresh) & (df.ttype == self.ttype)]
        return df[(df.signature == signature) & (df.treatment.isin(treatpool))]['treatment'].unique().tolist()

    def parse_exposure(self, signapool):

        expo = self.exposures[signapool]
        expo = expo[expo.index.isin(self.samples)]
        expo = expo.where(expo > 1, 0)  # small count drop to zero
        return expo

    def incidence(self, treatpool, signapool):
        """
        Returns:
            Boolean matrix with treatments and signatures matching treatments
            with signatures iff "treatment" is likely to elicit "signature"
        """
        df = self.discriminate
        df = df[(df.effect_size > self.effect_thresh) & (df.pvals < self.pvalue_thresh) & (df.ttype == self.ttype)]

        matrix = pd.DataFrame(index=treatpool, columns=signapool)
        for treat in treatpool:
            signatures = df[df['treatment'] == treat]['signature'].values
            for sig in signapool:
                matrix.loc[treat, sig] = (sig in signatures)
        return matrix

    def treatfit(self, treatpool, signapool):

        # exposures dataframe
        exposure = self.parse_exposure(signapool)

        # build treatment table
        treatment_table = self.treatments
        treatment_table = treatment_table[treatment_table.columns.intersection(list(treatpool))]
        treatment_table = treatment_table.loc[treatment_table.sum(axis=1) > 0, :]

        # select samples and make same samples in exposures and treatments
        samples = list(set(treatment_table.index).intersection(exposure.index))
        exposure = exposure.loc[samples, :]
        treatment_table = treatment_table.loc[samples, :]

        # NNLS: transform the input
        T = treatment_table.values
        E = exposure.values

        A, b = nnls_input(T, E)

        # NNLS: custom bounds -- some values must be zero
        lb = np.zeros(A.shape[1])  # lower bound
        incidence_matrix = self.incidence(treatpool, signapool)
        mask = np.reshape(incidence_matrix.values, (len(treatpool) * len(signapool), 1), order='F')
        mask = list(mask.flatten())
        eps = 1e-2
        ub = lb + np.inf
        ub[np.logical_not(mask)] = eps  # upper bound

        # NNLS: solver
        res = scipy.optimize.lsq_linear(A, b, bounds=(lb, ub))

        # format solution as a dataframe
        x_matrix = np.reshape(res.x, (T.shape[1], E.shape[1]), order='F')
        x_dataframe = pd.DataFrame(index=treatment_table.columns, columns=signapool)
        x_dataframe.loc[:, :] = x_matrix.astype(np.float32)

        # collapse near-zero coefficients
        x_dataframe = x_dataframe.where(x_dataframe > 1, 0)

        return x_dataframe


def origin(ttype, treatment, signature, exposures_path, discriminate_path, matrix_treatment_path):

    # determine efficiencies of competing treatments
    e = EfficiencyDrugSpecific(exposures_path, discriminate_path, matrix_treatment_path, ttype=ttype)

    # corner case 1
    if signature not in e.exposures.columns:
        return None, None

    # corner case 2
    if treatment not in e.treatments.columns:
        return None, None

    treatpool_list = e.select_pools(treatment, signature)

    # corner case 3
    if len(treatpool_list) == 0:
        return None, None

    fold_change_dict = {}  # we compute all the fold changes head-to-head
    for i, c in enumerate(treatpool_list):
        if c != treatment:
            x = e.treatfit([treatment] + [treatpool_list[i]], [signature])
            x = dict(zip(list(x.index), list(x.iloc[:, 0].values)))
            if x[c] == 0:
                fold_change_dict[c] = np.inf
            else:
                fold_change_dict[c] = x[treatment] / x[c]

    if len(fold_change_dict) == 0:
        etiology = treatment
        score = 1
    else:

        etiology = treatment
        for k, v in fold_change_dict.items():
            if v < 1:
                etiology = 'None'
        score = min(list(fold_change_dict.values()))
    return etiology, score


def dump_efficiency_dict(exposures_path, discriminate_path, matrix_treatment_path, output):

    eff = EfficiencyDrugSpecific(exposures_path, discriminate_path, matrix_treatment_path)
    with gzip.open(output, 'wb') as fn:
        pickle.dump(eff, fn)


def attribution(exposures_path, discriminate_path, matrix_treatment_path,  output):

    df_full = pd.read_csv(discriminate_path, sep='\t')
    toconcat = []
    for ttype, df in df_full.groupby(by='ttype'):
        df = df.astype({
            'signature': str, 'effect_size': float, 'pvals': float,
            'treatment': str, 'ttype': str, 'total_treated': float
        })

        df = df[(df.pvals < PVALUE) & (df.effect_size > EFFECT_SIZE)].copy()

        X = df.apply(lambda row: origin(
            row['ttype'], row['treatment'], row['signature'],
            exposures_path, discriminate_path, matrix_treatment_path
        ), axis=1)

        df['origin'] = list(map(lambda x: x[0], X.values))
        df['score'] = list(map(lambda x: x[1], X.values))

        toconcat.append(df)

    merged = pd.concat(toconcat)
    merged.to_csv(output, sep='\t', index=False)
