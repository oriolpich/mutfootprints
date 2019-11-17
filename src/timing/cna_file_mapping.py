# Import modules
from collections import defaultdict
from glob import glob
import gzip
import os

from tqdm import tqdm
import dill as pickle


def get_CNAS():
    # grep all cna files  and match it to the tumor ID
    all_cnas = glob('/workspace/datasets/hartwig/20190502/DR-024-update2/data/*/*purple.cnv')
    dic_cnas = defaultdict(str)

    for cna in tqdm(all_cnas):
        cna_s = os.path.basename(cna).split('.purple')[0]
        dic_cnas[cna_s] = cna

    os.makedirs("data/hartwig/clonal_structure/", exist_ok=True)
    pickle.dump(
        dict(dic_cnas),
        gzip.open('data/hartwig/clonal_structure/dic_cnas.pckl.gz', 'wb')
    )


if __name__ == '__main__':
    get_CNAS()
