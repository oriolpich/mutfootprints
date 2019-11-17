"""
Given a treatment, it gives a collection of signatures
that may be related to the treatment after adjusting by the tumor type.
"""

import gzip
import os
import sys

import numpy as np
import pandas as pd
import dill as pickle
from sklearn.linear_model import LogisticRegression
from scipy.stats import zscore

import warnings
warnings.filterwarnings('ignore')

scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)
from utils.metadata import return_metadata


np.random.seed(42)


def parse_treat(toregression, tumor_type, rand, all_sigs, raw_data):
    """
    Returns:
        filter is a list of tuples with colnames and treatment labels
    """

    treat_df = toregression.copy()

    if tumor_type != 'Pan':
        treat_df = treat_df[treat_df[tumor_type] == 1]

    treat_df = treat_df.transpose()
    treat_df = treat_df.loc[:, ~treat_df.columns.duplicated()]
    treat_df = treat_df.loc[treat_df.index.intersection(all_sigs + ['resp'])]  # index are signatures + treatment;
                                                                               # columns are samples;

    # filter low expression values
    totalexp = treat_df.sum(axis=1)
    totalexp = totalexp[totalexp > 10.]
    treat_df = treat_df.loc[list(totalexp.index) + ['resp'], :]

    # treatment/no-treatment is the response variable
    y = treat_df.loc['resp', :].values.astype(np.float64)
    # other treat_df values will be the design matrix
    treat_df.drop('resp', inplace=True)

    if not raw_data:

        # log + z-score transformation
        treat_df.loc[:, :] = np.log(treat_df.values + 1)
        treat_df.loc[:, :] = zscore(treat_df.values, axis=1)

    if filter is not None:
        cols, y = zip(*rand)
        y = np.array(list(map(lambda x: int(x.values[0]) if isinstance(x, pd.Series) else int(x), y)))
        treat_df = treat_df.loc[:, list(cols)]

    return treat_df, y


def run_logistic(to_regression, randomization, tumor_type, all_sigs, raw_data):

    treat_df, y = parse_treat(
        to_regression, tumor_type=tumor_type, rand=randomization, all_sigs=all_sigs, raw_data=raw_data
    )
    X = treat_df.values.astype(np.float64)
    logit = LogisticRegression(random_state=42, solver='lbfgs', penalty='l2', C=0.1).fit(X.T, y)
    dict_res = dict((zip(list(treat_df.index), logit.coef_[0])))
    return dict_res


class Urn:
    """
    Urn class provides with the container and methods to randomly sample
    balanced sets of treated vs untreated samples on a tumor-type by tumor-type basis.
    """

    def __init__(self, toregression, all_ttypes):

        dg = toregression[all_ttypes]
        self.encoding = {tt: tuple([int(tt == tt2) for tt2 in all_ttypes]) for tt in all_ttypes}
        dg['tumor_type'] = dg.apply(lambda r: self.get_tumor_type(r), axis=1)
        dg['resp'] = toregression["resp"].values
        self.bag = dg[['tumor_type', 'resp']]

    def get_tumor_type(self, row):

        row = tuple(row.values)
        for c in self.encoding:
            if row == self.encoding[c]:
                return c

    @staticmethod
    def balanced_sample(df):
        """sample a balanced set of indices with replacement"""

        df_one = df[df["resp"] == 1]
        df_zero = df[df["resp"] == 0]

        # split samples according to the number of treated patients.
        par = 3
        if len(df_one) < 50:

            par = 2

        m1, m2 = len(df_one), len(df_zero)
        n = min(m1, m2) // par

        if n > 0:
            ones = np.random.choice(list(df_one.index), size=n, replace=False)
            zeros = np.random.choice(list(df_zero.index), size=n, replace=False)
            samples = list(ones) + list(zeros)
        else:
            samples = []

        return [(s, df.loc[s, "resp"]) for s in samples]

    def stratified_random_sample(self, tumor_type='Pan'):
        """
        Generate a treated vs untreated balanced set of indices
        on a tumor-type by tumor-type basis
        """

        pool = []
        if tumor_type == 'Pan':
            for ttype in self.bag['tumor_type'].unique():
                pool += self.balanced_sample(self.bag[self.bag['tumor_type'] == ttype])
        else:
            return self.balanced_sample(self.bag[self.bag['tumor_type'] == tumor_type])
        return pool

    def randomize(self, n=1000, tumor_type='Pan'):

        np.random.seed(42)
        return [self.stratified_random_sample(tumor_type=tumor_type) for _ in range(n)]


def merge_dicts(dict_iterable):

    d = {}
    for k in dict_iterable[0].keys():
        d[k] = list(dict_[k] for dict_ in dict_iterable)
    return d


def effect_size(rand, to_regression, tumor_type, all_sigs):

    zeros = {}
    ones = {}
    effect_size_list = []

    for sample in rand:

        treat_df, y = parse_treat(to_regression, tumor_type=tumor_type, rand=sample, all_sigs=all_sigs, raw_data=True)
        df = treat_df.iloc[:, (y == 1)]
        for sign in df.index:
            ones[sign] = ones.get(sign, []) + df.loc[sign, :].values.tolist()

        dg = treat_df.iloc[:, (y == 0)]
        for sign in dg.index:
            zeros[sign] = zeros.get(sign, []) + dg.loc[sign, :].values.tolist()

        effect = {sign: fold_change(zeros[sign], ones[sign]) for sign in zeros}
        effect_size_list.append(effect)

    effect_dict = merge_dicts(effect_size_list)
    effect_dict = {k: np.mean(v) for k, v in effect_dict.items()}

    return effect_dict


def fold_change(zeros, ones):

    eps = 1
    m_zeros = np.mean(zeros)
    m_ones = np.mean(ones)  # rule out outlier values

    if m_zeros > eps:
        return m_ones / m_zeros
    else:
        return 0


def empirical_pvalue(values):

    assert (len(values) > 0)
    return len([v for v in values if v <= 0]) / len(values)


def create_matrix_treatments(drug, treated_path, exposures_path, tumor_type, keep_sigs=True):

    # load treated samples
    with gzip.open(treated_path, 'rb') as fd:
        treated_samples = pickle.load(fd)

    # return equivalence sample tumor type
    dic_primary_full, _ = return_metadata()

    # load mutational signature exposures at the Pan level
    d = pd.read_csv(exposures_path, sep='\t', index_col=0)

    # get all signatures
    all_sigs = d.index.tolist()

    # divide samples treated vs not treated
    if tumor_type == 'Pan':
        treated = list(treated_samples[drug]['Pan']['YES'])
        notreated = list(treated_samples[drug]['Pan']['NO'])
    else:
        treated = list(treated_samples[drug][tumor_type]['YES'])
        notreated = list(treated_samples[drug][tumor_type]['NO'])

    # select samples treated if they exist in the dataframe.
    s_treat = [s for s in d.columns if s in treated]
    s_nottreat = [s for s in d.columns if s in notreated]

    # get samples and transpose them
    affected = d[list(s_treat)].T
    notaffected = d[list(s_nottreat)].T

    affected_reg = affected.copy()
    notaffected_reg = notaffected.copy()

    # add variable response
    affected_reg['resp'] = 1
    notaffected_reg['resp'] = 0

    # concat all samples
    toregression = pd.concat([affected_reg, notaffected_reg])

    # if keep_sigs flag is True, we will select for the regression all the signatures which
    # have at least exposure in one fifth of the samples
    if keep_sigs:

        keep_exposed_signatures = []
        for signature in affected_reg.columns:
            exposed_to_sig = len(affected_reg[affected_reg[signature] > 10])
            if exposed_to_sig > (len(affected_reg) / 5):
                keep_exposed_signatures.append(signature)

        keep_exposed_signatures = keep_exposed_signatures + ['resp']
        toregression = toregression[keep_exposed_signatures]

    # add intercept
    toregression['_intercept'] = 1

    # add tumor type in the dataframe
    for ttype_l in set(dic_primary_full.values()):
        # get whether the sample belongs to a specific tumor type
        val = [1 if dic_primary_full.get(s, s) == ttype_l else 0 for s in toregression.index.tolist()]
        toregression[ttype_l] = val

    return toregression, all_sigs, list(set(list(dic_primary_full.values()))), len(treated)


def ensemble_regression(drug, treated_path, exposures_path,  ttype, n, keep_sigs=True):
    """
    This will run a regression to identify putative relations between drugs and treatments.
    It will generate n set of a balanced randomized splits.
    The size of the splits will depend on the num

    Args:
        drug: drug that we are analyzing
        treated_path: path where the treated/untreated dictionary is located
        exposures_path: path to signature extractions
        metadata_path: path to metadata
        output: path to output
        ttype: Tumor type we are analyzing. If None, it will do a Pan analysis
        n: Number of splits
        keep_sigs: Flag to remove uninformative signatures prior to the regression
    Returns:
        effect, pvals, merge_res, treated_samples: tuple with results of the regression analysis
    """

    # return the matrix, the signatures list and the tumor type list
    toregression, all_sigs, all_ttypes, treated_samples = create_matrix_treatments(
        drug, treated_path, exposures_path, ttype, keep_sigs
    )

    # randomize
    urn = Urn(toregression, all_ttypes)
    rand = urn.randomize(n=n, tumor_type=ttype)

    if len(rand[0]) == 0:
        raise ValueError('Randomization step cannot be done because of lack of treated samples')

    if len(all_sigs) < 2:
        raise ValueError('Randomization step cannot be done because of lack of putative signatures')

    # run logistic regression for each of the randomizations
    pool = []
    for randomization in rand:
        dict_res = run_logistic(toregression, randomization, tumor_type=ttype, all_sigs=all_sigs, raw_data=False)
        pool.append(dict_res)

    # get results
    merge_res = merge_dicts(pool)

    # get the effect size
    effect = effect_size(rand, toregression, ttype, all_sigs)  # fold change
    pvals = {k: empirical_pvalue(v) for k, v in merge_res.items()}

    return effect, pvals, merge_res, treated_samples


