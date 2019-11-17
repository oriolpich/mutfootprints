# this is to decompose some drug names or combinations
composed_drugs = {

    'Docetaxel, Doxorubicin, Cyclophosphamide (TAC)': [
        'Docetaxel', 'Doxorubicin', 'Cyclophosphamide'
    ],
    'Dexamethasone, High-dose cytarabine, Cisplatin (DAHP) - Ifosfamide, mitoxantrone and etoposide (VIM) - Dexamethasone, High-dose cytarabine, Cisplatin (DAHP)': [
        'Dexamethasone', 'Cytarabine', 'Cisplatin', 'Ifosfamide', 'Mitoxantrone', 'Etoposide'
    ],
    'Fluorouracil, Epirubicin, Cyclophosphamide, Docetaxel (FEC-D)': [
        'Fluorouracil', 'Epirubicine', 'Cyclophosphamide', 'Docetaxel'
    ],
    'Fluorouracil, Epirubicine and Cyclofosfamide (FEC)': [
        'Fluorouracil', 'Epirubicine', 'Cyclofosfamide'
    ],
    'Fluorouracil, Leucovorin Calcium (FU': [
        'Fluorouracil', 'Leucovorin'
    ],
    'Tegafur, Gimeracil, Oteracil (S1)': [
        'Tegafur', 'Gimeracil', 'Oteracil'
    ],
    'Procarbazine, Lomustine, Vincristine (PCV)': [
        'Procarbazine', 'Lomustine', 'Vincristine'
    ],
    'Doxorubicin, Bleomycin, Vinblastine, Dacarbazin (ABVD)': [
        'Doxorubicin', 'Bleomycin', 'Vinblastine', 'Dacarbazin'
    ],
    'Doxorubicin, Cyclophosphamide (AC)': [
        'Doxorubicin', 'Cyclophosphamide'
    ],
    'Pemetrexed, Carboplatin': [
        'Premetrexed', 'Carboplatin'
    ],
    'leucovorin, fluorouracil, irinotecan (FOLFIRI)': [
        'Leucovorin', 'Fluorouracil', 'Irinotecan'
    ],
    'leucovorin, fluorouracil, oxaliplatin (FOLFOX)': [
        'Leucovorin', 'Fluorouracil', 'Oxaliplatin'
    ],
    'mustine, Vincristine Sulfate, Procarbazine Hydrochloride, Prednisone (MOPP)': [
        'Limustine', 'Vincristine', 'Procarbazine', 'Prednisone'
    ],
    'Vincristine Sulfate, Dactinomycin, Cyclophosphamide (VAC)': [
        'Vincristine', 'Dactinomycin', 'Cyclophosphamide'
    ],
    'Vincristine Sulfate, Etoposide Phosphate, Prednisone, Doxorubicin Hydrochloride (OEPA)': [
        'Vincristine', 'Etoposide', 'Prednisone', 'Doxorubicin'
    ],
    'Cabazitaxel, Zoledronic Acid': [
        'Cabazitaxel', 'Zoledronic Acid'
    ],
    'Capecitabine, Oxaliplatin (CAPOX)': [
        'Capecitabine', 'Oxaliplatin'
    ],
    'Carboplatin, Paclitaxel (PC)': [
        'Carboplatin', 'Paclitaxel'
    ],
    'Carmustine, Etoposide, Cytarabine, Melphalan (BEAM)': [
        'Carmustine', 'Etoposide', 'Cytarabine', 'Melphalan'
    ],
    'Cyclophosphamide, Methotrexate, Fluorouracil (CMF)': [
        'Cyclophosphamide', 'Methotrexate', 'Fluorouracil'
    ],
    'Cyclophosphamide, Vincristine Sulfate, Procarbazine Hydrochloride, Prednisone (COPP)': [
        'Cyclophosphamide', 'Vincristine', 'Procarbazine', 'Prednisone',
    ],
    'Rituximab, Cyclophosphamide, Doxorubicin, Vincristine, Prednisolone (R-CHOP)': [
        'Rituximab', 'Cyclophosphamide', 'Doxorubicin', 'Vincristine', 'Prednisolone'
    ],
    'ABVD': [
        'Doxorubicin', 'Bleomycine', 'Vinblastine', 'Dacarbazine'
    ],
    'AC': [
        'Doxorubicin', 'Cyclophosphamide'
    ],
    'Ad[1': [
        'Doxorubicin'
    ],
    'Adriamycin': [
        'Doxorubicin'
    ],
    'BEAM': [
        'Carmustine', 'Etoposide', 'Cytarabine', 'Melphalan'
    ],
    'CMF': [
        'Cyclophosphamide', 'Methotrexate', 'Fluorouracil'
    ],
    'COPP': [
        'Cyclophosphamide', 'Vincristine', 'Procarbazine', 'Prednisone'
    ],
    'Capox': [
        'Capecitabine', 'Oxaliplatin'
    ],
    'DHAP-VIM-DHAP': [
        'Rituximab', 'Dexamethason', 'Cytarabine', 'Cisplatin', 'Etoposide', 'Ifosfamide', 'Methotrexate'
    ],
    'FEC': [
        'Fluorouracil', 'Epirubicine', 'Cyclophosphamide'
    ],
    'FEC-D': [
        'Fluorouracil', 'Epirubicine', 'Cyclophosphamide', 'Docetaxel'
    ],
    'Folfirinox': [
        'Fluorouracil', 'Irinotecan', 'Oxaliplatin'
    ],
    'MOPP': [
        'Nitrogen_mustard', 'Oncovin', 'Procarbazine', 'Prednisone'
    ],
    'OEPA': [
        'Oncovin', 'Etoposide', 'Prednisone', 'Doxorubicin'
    ],
    'PCV': [
        'Procarbazine', 'Oncovin', 'Lomustine'
    ],
    'R-CHOP': [
        'Rituximab', 'Cyclophosphamide', 'Doxorubicin', 'Oncovin', 'Prednisolone'
    ],
    'S1': [
        'Fluorouracil'
    ],  # plus others that enhance its activity
    'SYD985': [
        'Trastuzumab'
    ],
    'TAC': [
        'Docetaxel', 'Doxorubicin', 'Cyclophosphamide'
    ],
    'Taxol': [
        'Paclitaxel'
    ],
    'VAC': [
        'Oncovin', 'Dactinomycin', 'Cyclophosphamide'
    ],
    'Zoladex': [
        'Gosereline'
    ],
    'Casodex': [
        'Bicalutamide'
    ],
    'Eligard': [
        'Leuproreline'
    ],
    'Carbotaxol': [
        'Paclitaxel', 'Carboplatin'
    ],
    'LV': [
        'Fluorouracil', 'Folic acid'
    ],
    'Abirateron, Prednison': [
        'Abirateron', 'Prednisone'
    ],
    'Carboplatin, Gemcitabine': [
        'Carboplatin', 'Gemcitabine'
    ],
    'GDC': [
        'Gemcitabine', 'Docetaxel', 'Carboplatin'
    ],
    'GDC-0032': [
        'Gemcitabine', 'Docetaxel', 'Carboplatin'
    ],
    'Letrozol (+ LEE01': [
        'Letrozole', 'LEE01'],
    'Fluorouracil, Epirubicin, Cyclophosphamide (FEC)': [
        'Fluorouracil', 'Epirubicin', 'Cyclophosphamide'
    ],
    'Mechlorethamine Hydrochloride, Vincristine Sulfate, Procarbazine Hydrochloride, Prednisone (MOPP)': [
        'Mechlorethamine', 'Vincristine', 'Procarbazine', 'Prednisone'
     ],
    'Doxorubicin, Bleomycin, Vinblastine, Dacarbazine (ABVD)': [
        'Doxorubicin', 'Bleomycin', 'Vinblastine', 'Dacarbazine'
    ],
    'Leucovorin, fluorouracil, irinotecan (FOLFIRI)': [
        'Leucovorin', 'Fluorouracil', 'Irinotecan'
    ],
    'Abiraterone, Prednisone': [
        'Abiraterone', 'Prednisone'
    ],
    'Letrozole (+ LEE011 or placebo in Monaleesa trial)': [
      'Letrozole', 'Unknown'
    ],
    'Leucovorin, fluorouracil, oxaliplatin (FOLFOX)': [
        'Leucovorin', 'Fluorouracil', 'Oxaliplatin'
    ]

}

# add instances not found in the FDA dict
drugs_class_missing = {

    'UNK': 'Unknown',
    'UNKNOWN CHEMOTHERAPY': 'Unknown',
    'PLACEBO': 'Placebo',
    'RIBOCICLIB OR PLACEBO': 'Unknown',
    'PEMBROLIZUMAB OR PLACEBO': 'Unknown',
    'TUMOR TREATING FIELDS (TTFIELDS)': 'Experimental therapy',
    'UNKNOWN': 'Unknown',
    'COMPLETE HEMATOLOGIC RESPONSE (CHR) 3996 TRIAL': 'Experimental therapy',
    'PLACEBO IN MONALEESA TRIAL)': 'Placebo',
    'STUDY DRUG (RIBOCILIB OR PLACEBO)': 'Unknown',
    'ACTINIUM-225 PROSTATE-SPECIFIC MEMBRANE ANTIGEN (PSMA) RADIOLIGANDTHERAPIE': 'Nuclear therapy',
    'ABIRATERON': 'Androgen Receptor Inhibitor',
    'ETHINYLESTRADIOL': 'Estrogen Receptor Agonist',
    'ADT': 'Androgen Receptor Inhibitor',
    'AMG 211': 'Autologous Cellular Immunotherapy',
    'AMG232': 'Miscellanious',
    'ACTINIUM-225 PSMA RADIOLIGANDTHERAPIE': 'Nuclear therapy',
    'ANDROCUR': 'Androgen Receptor Inhibitor',
    'ANTHRACYCLINES': 'Anthracycline Topoisomerase Inhibitor',
    'ANTHRACYCLINE': 'Anthracycline Topoisomerase Inhibitor',
    'ANTI-GITR': 'Autologous Cellular Immunotherapy',
    'ANTI - GLUCOCORTICOID INDUCED TNF RECEPTOR (GITR)': 'Autologous Cellular Immunotherapy',
    'BOLLA HORMONE TREATMENT': 'Miscellanious',
    'AROMATASE INHIBITOR': 'Aromatase Inhibitor',
    'AROMATASE INHIBITOR (AI)': 'Aromatase Inhibitor',
    'BEZ235': 'Kinase inhibitor',
    'ALPELISIB': 'Kinase inhibitor',
    'INC280': 'Kinase Inhibitor',
    'BMS 986148': 'Mesothelin',
    'BMS986156': 'Autologous Cellular Immunotherapy',
    'BINIMETINIB': 'Kinase Inhibitor',
    'BRACHYTHERAPY': 'Radiotherapy',
    'CHR-3996 TRIAL': 'Histone Deacetylase Inhibitor',
    'CABAZITAXEL': 'Microtubule Inhibitor',
    'CYPROTERON': 'Androgen Receptor Inhibitor',
    'CYPROTERONE': 'Androgen Receptor Inhibitor',
    'DMOT4039A': 'Mesothelin',
    'DEGARELIX': 'Gonadotropin Releasing Hormone Antagonist',
    'DENDRITIC CELL THERAPY': 'Autologous Cellular Immunotherapy',
    'DENDRITIC CELL THERAPY (DCVAC)': 'Autologous Cellular Immunotherapy',
    'ENCORAFENIB': 'Kinase Inhibitor',
    'FOLIC ACID': 'Folinic Acid',
    'HORMONAL THERAPY': 'Aromatase Inhibitor',
    'HYPERTHERMIA': 'Hyperthermia',
    'IGF INHIBITOR (STUDY)': 'Miscellanious',
    'INC 280': 'Kinase Inhibitor',
    'LHRH AGONIST': 'Gonadotropin Releasing Hormone Receptor Agonist',
    'MYELOID DENDRITICAL CELLS': 'Autologous Cellular Immunotherapy',
    'OLARATUMAB': 'PDGF antibody',
    'ORTERONEL': 'Aromatase Inhibitor',
    'ORTERONEL ': 'Aromatase Inhibitor',
    'PF-04518600': 'Autologous Cellular Immunotherapy',
    'PPT-E1A]': 'Oncolytic adenovirus',
    'RFA': 'Radiofrequency',
    'RFASE': 'Radiofrequency',
    'RADIUM-223': 'Nuclear therapy',
    'RILIMOGENE': 'Autologous Cellular Immunotherapy',
    'ROCILETINIB': 'Kinase Inhibitor',
    'SAMARIUM-153': 'Nuclear therapy',
    'SUPREFACT': 'Gonadotropin Releasing Hormone Receptor Agonist',
    'T-VEC (TALIMOGEEN)': 'Oncolytic adenovirus',
    'TTFIELDS': 'Experimental therapy',
    'TASELISIB': 'Kinase Inhibitor',
    'THALIDOMIDE': 'Thalidomide Analog',
    'TRANSARTERIAL CHEMO-EMBOLIZATION': 'Unknown',
    'UNKNOWN (NEOADJUVANT CHEMOTHERAPY)': 'Unknown',
    'UNKNOWN ADJUVANT CHEMOTHERAPY': 'Unknown',
    'UNKNOWN CHEMOTHERAPY FOR BREAST TUMOR': 'Unknown',
    'VELIPARIB': 'Poly(ADP-Ribose) Polymerase Inhibitor',
    'HVEGF': 'Vascular Endothelial Growth Factor Receptor 2 Antagonist',
    'NAN': 'Unknown',
    'PLACEBO OR TASELISIB': 'Unknown',
    'PIMASERTIB': 'Kinase Inhibitor',
    ' PLACEBO (PT HAD PLACEBO)': 'PLACEBO',
    'APALUTAMIDE': 'Androgen Receptor Inhibitor',
    'BUSERELIN': 'Gonadotropin Releasing Hormone Receptor Agonist',
    'ANDROGEN DEPRIVATION THERAPY (ADT)': 'Androgen Receptor Inhibitor',
    'ANTI-ANDROGEN': 'Androgen Receptor Inhibitor',
    'ANTIANDROGEN': 'Androgen Receptor Inhibitor',
    'LUTEINIZING HORMONE-RELEASING HORMONE (LHRH) AGONISTS': 'Gonadotropin Releasing Hormone Receptor Agonist',
    'COMPLETE HEMATOLOGIC RESPONE (CHR) 3996 TRIAL': 'Miscellanious',
    'DENDRITIC CELL VACCINATION (DCVAC) OR PLACEBO': 'Miscellanious',
    'DENDRITIC CELLS OR PLACEBO (STUDY)': 'Miscellanious',
    'MYELOID DENDRITIC CELLS': 'Autologous Cellular Immunotherapy',
    'TREMELIMUMAB': 'CTLA-4-directed Blocking Antibody',
    'NIVOLUMAB OR PLACEBO': 'Unknown',
    'DENDRITIC CELL ACTIVATION OR PLACEBO': 'Miscellanious',
    'LEE011': 'Kinase Inhibitor',
    'PLACEBO VS OLAPARIB': 'Placebo',
    'LEVEREMBOLISATIE': 'Unknown',
    'HYPERTHERMIC INTRAPERITONEAL CHEMOTHERAPY (HIPEC) - COLOPEC': 'Unknown',
    'TALIMOGENE LAHERPAREPVEC (T-VEC)': 'Oncolytic adenovirus',
    'OTERACIL': 'Unknown',
    'GIMERACIL': 'Unknown',
    'PLACEBO OR RIBOCICLIB': 'Unknown',
    'TAS-119': 'Kinase Inhibitor',
    'DENDTRITIC CELACTIVATION OR PLACEBO': 'Unknown',
    'LEE0011 OR PLACEBO': 'Unknown',
    'HYPERTHERMIC INTRAPERITONEAL CHEMOTHERAPY (HIPEC)': 'Unknown',
    'UNKNOW': 'Unknown',
    'RADIUM OR PLACEBO': 'Unknown',
    'INSULIN-LIKE GROWTH FACTOR (IGF) INHIBITOR (STUDY)': 'Insulin-like growth factor inhibitor',
    'LEE01': 'Kinase Inhibitor',
    'PREMETREXED': 'Folate antagonist',
    'TUMOR NECROSIS FACTOR ALPHA': 'Unknown',
    'VINDESINE': 'Vinca Alkaloid',
    'LUCANIX': 'Autologous Cellular Immunotherapy',
    'RADIUS INVESTIGATIONAL DRUG ELACESTRANT (RAD1901)': 'Estrogen Receptor Antagonist',
    'SOMATOSTATINE': 'Somatostatin Analog',
    'TIPIRACIL': 'Nucleoside Metabolic Inhibitor',
    'RADIOFREQUENCY ABLATION (RFA)': 'Radiofrequency',
    'SHP-2 INHIBITOR': 'Unknown',
    'SSA': 'Unknown',
    'SAR566658': 'Autologous Cellular Immunotherapy',
    'MITOXANTRONE': 'Topoisomerase Inhibitor',
    'LORLATINIB': 'ROS1 Inhibitor',
    'PEPTIDE RECEPTOR RADIONUCLIDE THERAPY (PRRT)': 'Nuclear therapy',
    'LURBINECTEDIN': 'Alkylating Drug',
    'ISOLATED LIMB PERFUSION (ILP) RIGHT LEG': 'Unknown',
    'VASCULAR ENDOTHELIAL GROWTH FACTOR (VEGF) - INHIBITOR': 'Vascular Endothelial Growth Factor Receptor 2 Antagonist',
    'COBIMETINIB': 'Kinase Inhibitor',
    'PLACEBO OR OLAPARIB': 'Unknown',
    'ORTERONEL OR PLACEBO': 'Unknown',
    'DENDRITIC CELL THERAPY (DCVAC) OR PLACEBO': 'Unknown',
    'ATEZOLIZUMAB OR PLACEBO': 'Unknown',
    'AZD3514': 'Androgen Receptor Inhibitor',
    'LEE011 OR PLACEBO': 'Unknown',
    'BMS 986156': 'Autologous Cellular Immunotherapy',
    'LIVER EMBOLIZATION': 'Unknown',
    'GDC OR PLACEBO': 'Unknown',
    'ACTINIUM-225 PROSTATE-SPECIFIC MEMBRANE ANTIGEN (PSMA) RADIOLIGAND THERAPY': 'Nuclear therapy',
    'GDC-0032 OR PLACEBO': 'Unknown',
    'OLARATUMAB OR PLACEBO': 'Unknown',
    'RADIUM-223 OR PLACEBO': 'Unknown',
    'TRANSARTERIAL CHEMOEMBOLIZATION': 'Unknown',
    'AMG 232': 'Experimental therapy',
    'LUTETIUM': 'Nuclear therapy',
    'IODINE-131': 'Nuclear therapy',
    'CORTICOSTEROID': 'Corticosteroid',
    'GALLIUM': 'Nuclear therapy',

}

# fix some of the input names
fix_input_names = {
    'Doxorubicine': 'DOXORUBICIN HYDROCHLORIDE',
    'Doxorubicin': 'DOXORUBICIN HYDROCHLORIDE',
    'Bleomycine': 'BLEOMYCIN SULFATE',
    'Bleomycin': 'BLEOMYCIN SULFATE',
    'Doxorubicin (liposomal)': 'DOXORUBICIN HYDROCHLORIDE',
    'Dexamethason': 'DEXAMETHASONE SODIUM PHOSPHATE',
    'Dexamethasone': 'DEXAMETHASONE SODIUM PHOSPHATE',
    'Orteronel ': 'ORTERONEL',
    'Gosereline': 'GOSERELIN ACETATE',
    'Irinotecan': 'IRINOTECAN HYDROCHLORIDE',
    'Leuproreline': 'LEUPROLIDE ACETATE',
    'Leuprorelin': 'LEUPROLIDE ACETATE',
    'Nitrogen_mustard': 'MECHLORETHAMINE HYDROCHLORIDE',
    'Oncovin': 'VINCRISTINE SULFATE',
    'Procarbazine': 'PROCARBAZINE HYDROCHLORIDE',
    'Vinblastine': 'VINBLASTINE SULFATE',
    'Vincristine': 'VINCRISTINE SULFATE',
    'Abiraterone': 'ABIRATERONE ACETATE',
    'BCG': 'BACILLUS CALMETTE-GUERIN SUBSTRAIN TICE LIVE ANTIGEN',
    'Bacillus Calmette-Guerin (BCG)': 'BACILLUS CALMETTE-GUERIN SUBSTRAIN TICE LIVE ANTIGEN',
    'Bacillus Calmette-Guérin (BCG)': 'BACILLUS CALMETTE-GUERIN SUBSTRAIN TICE LIVE ANTIGEN',
    'Bisphosphonates': 'ALENDRONATE SODIUM',
    'Bisphosphonate': 'ALENDRONATE SODIUM',
    'Busulphan': 'BUSULFAN',
    'Cabozantinib': 'CABOZANTINIB S-MALATE',
    'Chloroquine': 'HYDROXYCHLOROQUINE SULFATE',
    'Dabrafenib': 'DABRAFENIB MESYLATE',
    'Epi-adriamycine': 'DOXORUBICIN HYDROCHLORIDE',
    'Eribuline': 'ERIBULIN MESYLATE',
    'Eribulin': 'ERIBULIN MESYLATE',
    'Estramustine': 'ESTRAMUSTINE PHOSPHATE SODIUM',
    'FU': 'FLUOROURACIL',
    'Fluoropyrimidine': 'FLUOROURACIL',
    'Folinic acid': 'FOLIC ACID',
    'HIPEC': 'Unknown',
    'HIPEC-COLOPEC': 'Unknown',
    'Imatinib': 'IMATINIB MESYLATE',
    'Lapatinib': 'LAPATINIB DITOSYLATE',
    'Liposomal doxorubicin': 'DOXORUBICIN HYDROCHLORIDE',
    'Lucrin': 'LEUPROLIDE ACETATE',
    'Mitomycin-C': 'MITOMYCIN',
    'Mitomycin C': 'MITOMYCIN',
    'Pegylated liposomal doxorubicin': 'DOXORUBICIN HYDROCHLORIDE',
    'Pemetrexed': 'PEMETREXED DISODIUM HEMIPENTAHYDRATE',
    'Regorafenib': 'REGORAFENIB MONOHYDRATE',
    'Sandostatine': 'OCTREOTIDE ACETATE',
    'Somatuline': 'LANREOTIDE ACETATE',
    'Sunitinib': 'SUNITINIB MALATE',
    'TDM-1': 'TRASTUZUMAB',
    'Tamoxifen': 'TAMOXIFEN CITRATE',
    'Trametinib': 'TRAMETINIB DIMETHYL SULFOXIDE',
    'Triptorelin': 'TRIPTORELIN PAMOATE',
    'Ado-trastuzumab emtansine (TDM-1)': 'TRASTUZUMAB',
    'Anastrozol': 'ANASTROZOLE',
    'Chloromethine': 'ESTRAMUSTINE PHOSPHATE SODIUM',
    'Cytarabin': 'CYTARABINE',
    'Dacarbazin': 'DACARBAZINE',
    'Dexamethason (liposomal)': 'DEXAMETHASONE SODIUM PHOSPHATE',
    'Dexamethasone (liposomal)': 'DEXAMETHASONE SODIUM PHOSPHATE',
    'Leucovorin': 'LEUCOVORIN CALCIUM',
    'Tegafur': 'FLUOROURACIL',
    'Gemcitabine': 'GEMCITABINE HYDROCHLORIDE',
    'Cyclofosfamide': 'CYCLOPHOSPHAMIDE',
    'Docetaxel (CriPec)': 'DOCETAXEL',
    'Letrozol': 'LETROZOLE',
    'Interferon': 'INTERFERON GAMMA-1B',
    'Medroxyprogesteron': 'MEDROXYPROGESTERONE ACETATE',
    'Vinorelbin': 'VINORELBINE TARTRATE',
    'LV)': 'LEUCOVORIN CALCIUM',
    'Limustine': 'LOMUSTINE',
    'Nab-paclitaxel': 'PACLITAXEL',
    'Epirubicine': 'EPIRUBICIN HYDROCHLORIDE',
    'Epirubicin': 'EPIRUBICIN HYDROCHLORIDE',
    'Ado-trastuzumab emtansine (T-DM1)': 'TRASTUZUMAB',
    'Somatostatin': 'SOMATOSTATINE',
    'LU-177-PSMA': 'LUTETIUM',
    'LUTETIUM-177 DOTATATE': 'LUTETIUM',
    'Lutetium (Lu) -177- Prostate-specific membrane antigen (PSMA)': 'LUTETIUM',
    'iodine (I) -131': 'IODINE-131',
    'Lutetium (Lu) -177 DOTATATE':  'LUTETIUM',
    'Gallium(Ga) 68-DOTATATE': 'GALLIUM',
    'GALLIUM-68 (68GA)-DOTATATE': 'GALLIUM',
    'I-131 5500 MBQ': 'IODINE-131',
    'RADIOACTIVE IODINE': 'IODINE-131',
    '68GA-DOTATATE': 'GALLIUM',
    'Prednisone': 'CORTICOSTEROID',
    'RADIUM-223 DICHLORIDE': 'Nuclear therapy',
}
