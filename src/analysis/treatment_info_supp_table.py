import pandas as pd
import gzip
import pickle
from collections import defaultdict

FDA_treats = pickle.load(gzip.open('data/clinical_data/hartwig_FDA_treatments.pckl.gz'))
dic_treated = defaultdict(dict)
for drug, d in FDA_treats.items():
    if drug not in ['RADIATION', 'TOPOII']:
        for ttype, di in d.items():
            if ttype !='Pan':
                treated_s = len(di['YES'])
                dic_treated[drug][ttype] = treated_s

pd.DataFrame(dic_treated).fillna(0).to_csv('figures/table_patients_treated.tsv', sep ='\t',
                                           index = True, header = True)