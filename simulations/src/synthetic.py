import os
import json
import pickle
import gzip

import numpy as np
import scipy
import pandas as pd
from scipy.stats import beta
import matplotlib.pyplot as plt


# config supplementary table

supp_table_path = os.path.join('../data', 'contrib_mut_burden_signature_analyzer.csv')


class Signatures:

    def __init__(self, signatures_path):

        with gzip.open(signatures_path, 'rb') as f:
            self.table = pickle.load(f)
        self.channel_index = np.array(self.table.index)
        self.etiologies = {'platinum': 'SBS31', 'poleta': 'SBS9'}


def only_platinum(r):

    if (r['CARBOPLATIN'] + r['CISPLATIN'] + r['OXALIPLATIN'] > 0) and (r['CAPECITABINE/5-FU'] == 0):
        return True
    return False


def total(r):

    try:
        return r['treatment'] / r['proportion']
    except ZeroDivisionError:
        print(r['treatment'], r['proportion'])


def proportion_platinum_mutations():

    # load supplementary table

    df = pd.read_csv(supp_table_path, sep='\t')
    df.fillna(0, inplace=True)
    df.rename(columns={'Proportion treatment related mutations': 'proportion'}, inplace=True)
    df.rename(columns={'Total treatment related mutations': 'treatment'}, inplace=True)

    # discard high total burden >= 1e5

    df = df[df['proportion'] > 0]
    df['only_platinum'] = df.apply(only_platinum, axis=1)
    df_platinum = df[df['only_platinum']].copy()
    df_platinum['total'] = df_platinum.apply(total, axis=1)
    df_platinum = df_platinum.loc[lambda x: x.total < 1e5]

    # return proportion values

    prop_platinum = df_platinum[['proportion']].copy()
    return prop_platinum.values[:, 0]


class Synthetic:

    def __init__(self, act, proc):

        self.sample_index  = np.array(act.columns)
        self.process_index = np.array(act.index)
        self.channel_index = np.array(proc.index)
        self.act  = act.values                    # cols=samples; rows=processes
        self.proc = proc[list(act.index)].values  # cols=processes; rows=channels
        self.params = beta.fit(proportion_platinum_mutations())
        self.etiologies = {'platinum': 'SBS31', 'poleta': 'SBS9'}

    @property
    def templates(self):

        channel_exposure = np.dot(self.proc, self.act)

        return dict(zip(self.sample_index, list(channel_exposure.T)))

    def pick_samples(self, n_samples):

        return list(np.random.choice(self.sample_index, size=n_samples, replace=True))

    def primary_draw(self, n_samples):

        df_dict = {}
        samples = self.pick_samples(n_samples)
        for i, s in enumerate(samples):
            template = self.templates[s]
            total = np.sum(template)
            profile = template / total
            df_dict[f'sample_{i+1}_{s.split("__")[0]}'] = list(np.random.multinomial(int(total), pvals=profile))

        return pd.DataFrame(df_dict)

    def proportion_treated(self, n):

        """
        self.params were obtained upon fitting the estimated
        proportions of mutations attributable to platinum-based drugs
        according to table in global path supp_table_path
        """

        a, b, loc, scale = self.params

        return beta.rvs(a, b, loc=loc, scale=scale, size=n)

    def full_draw(self, n_samples, n_treated, sig):

        catalogue = self.primary_draw(n_samples)
        catalogue.rename(columns=dict(zip(list(catalogue.columns)[:n_treated],
                                          list(map(lambda s: s + '_t',
                                                   list(catalogue.columns)[:n_treated])))),
                         inplace=True)

        beta_p = self.proportion_treated(n_samples)
        all_burden = catalogue.sum(axis=0).values
        burden = np.array([b*(beta_p[i]/(1 - beta_p[i])) if i < n_treated else 0 for i, b in enumerate(all_burden)])
        index_sig = list(self.process_index).index(sig)
        profile = self.proc[:, index_sig]

        df_dict = {}
        for i, s in enumerate(catalogue.columns):
            df_dict[s] = list(np.random.multinomial(int(burden[i]), pvals=profile))
        treat_df = pd.DataFrame(df_dict)
        catalogue.loc[:, :] += treat_df.loc[:, :]
        res = {'catalogue': catalogue, 'burden': burden}

        return res


def generate_synthetic_data(synthetic, output_folder):

    np.random.seed(42)
    for sig in ['SBS9', 'SBS31']:
        for n_treated in [10, 25, 50, 100, 150]:
            for replicate in range(25):
                res = synthetic.full_draw(300, n_treated, sig)
                res['catalogue'].to_csv(os.path.join(output_folder,
                                                     f'{sig}.{n_treated}.{replicate}.catalogue.tsv'),
                                        sep='\t', index=False)
                with open(os.path.join(output_folder,
                                       f'{sig}.{n_treated}.{replicate + 1}.burden.json'), 'wt') as f:
                    json.dump(list(res['burden']), f)


def test_beta_fit():

    a, b, loc, scale = beta.fit(proportion_platinum_mutations())
    rvs = beta.rvs(a, b, loc=loc, scale=scale, size=1000)
    plt.hist(rvs, bins=50)
    plt.show()


if __name__ == '__main__':

    """
    python synthetic.py
    """

    test_beta_fit()
