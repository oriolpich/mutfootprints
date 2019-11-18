import os
import json
import numpy as np
import pandas as pd
from sklearn import metrics


def signature_recovery(sig, n_treated, replicate, signatures, deconstruction_folder):

    twin_signals = False
    case = f'{sig}.{n_treated}.{replicate}'
    df = pd.read_csv(os.path.join(deconstruction_folder,
                                  f'{case}/{case}.processes.tsv'), sep='\t')
    cosines = []
    for col in df.columns:
        reference = signatures[sig].values
        profile = df[col].values
        cosine = metrics.pairwise.cosine_similarity([reference], [profile])[0][0]
        cosines.append(cosine)
    similar = [c for c in cosines if c > 0.8]
    if len(similar) > 1:
        twin_signals = True
    cosine = max(cosines)
    index = [np.argmax(cosines)]
    if len(similar) > 1:
        index += list(np.argsort(cosines)[::-1][1: len(similar)])
    return index, cosine, twin_signals


def exposure_recovery(sig, n_treated, replicate, signatures, catalogue_folder, deconstruction_folder):

    # load synthetic burden

    burden_fn = f'{sig}.{n_treated}.{replicate + 1}.burden.json'
    with open(os.path.join(catalogue_folder, burden_fn), 'rt') as f:
        burden = json.load(f)
    burden = np.array(burden[:n_treated])

    # load exposures

    case = f'{sig}.{n_treated}.{replicate}'
    exposure = pd.read_csv(os.path.join(deconstruction_folder, f'{case}', f'{case}.exposures.tsv'), sep='\t')
    index, cosine, twin_signals = signature_recovery(sig, n_treated, replicate, signatures, deconstruction_folder)
    exposure = exposure.iloc[index[0], :n_treated].values

    return burden, exposure


def exposure_recovery_nontreated(sig, n_treated, replicate, signatures, catalogue_folder, deconstruction_folder):

    # load exposures

    case = f'{sig}.{n_treated}.{replicate}'
    exposure = pd.read_csv(os.path.join(deconstruction_folder, f'{case}/{case}.exposures.tsv'), sep='\t')
    index, cosine, twin_signals = signature_recovery(sig, n_treated, replicate, signatures, deconstruction_folder)
    exposure = exposure.iloc[index[0], n_treated + 1:].values

    return exposure


def concordance_cc(y_true, y_pred):

    """
    Read more: https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
    Original paper: Lawrence, I., and Kuei Lin. "A concordance correlation coefficient
    to evaluate reproducibility." Biometrics (1989): 255-268.
    Credit: https://github.com/stylianos-kampakis/supervisedPCA-Python/blob/master/Untitled.py

    Args:
        y_true : array-like. Ground truth (correct) target values.
        y_pred : array-like. Estimated target values.
    Returns:
       loss : float in the range [-1,1].
    """

    cor = np.corrcoef(y_true, y_pred)[0][1]

    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)

    var_true = np.var(y_true)
    var_pred = np.var(y_pred)

    sd_true = np.std(y_true)
    sd_pred = np.std(y_pred)

    numerator = 2 * cor * sd_true * sd_pred

    denominator = var_true + var_pred + (mean_true - mean_pred) ** 2

    return numerator / denominator
