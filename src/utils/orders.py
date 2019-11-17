
import numpy as np

def dbs_types():

    types = [
        'AC', 'AC', 'AC', 'AC', 'AC', 'AC', 'AC', 'AC', 'AC', 'AT', 'AT',
        'AT', 'AT', 'AT', 'AT', 'CC', 'CC', 'CC', 'CC', 'CC', 'CC', 'CC',
        'CC', 'CC', 'CG', 'CG', 'CG', 'CG', 'CG', 'CG', 'CT', 'CT', 'CT',
        'CT', 'CT', 'CT', 'CT', 'CT', 'CT', 'GC', 'GC', 'GC', 'GC', 'GC',
        'GC', 'TA', 'TA', 'TA', 'TA', 'TA', 'TA', 'TC', 'TC', 'TC', 'TC',
        'TC', 'TC', 'TC', 'TC', 'TC', 'TG', 'TG', 'TG', 'TG', 'TG', 'TG',
        'TG', 'TG', 'TG', 'TT', 'TT', 'TT', 'TT', 'TT', 'TT', 'TT', 'TT',
        'TT'
    ]

    subtypes = [
        'CA', 'CG', 'CT', 'GA', 'GG', 'GT', 'TA', 'TG', 'TT', 'CA', 'CC',
        'CG', 'GA', 'GC', 'TA', 'AA', 'AG', 'AT', 'GA', 'GG', 'GT', 'TA',
        'TG', 'TT', 'AT', 'GC', 'GT', 'TA', 'TC', 'TT', 'AA', 'AC', 'AG',
        'GA', 'GC', 'GG', 'TA', 'TC', 'TG', 'AA', 'AG', 'AT', 'CA', 'CG',
        'TA', 'AT', 'CG', 'CT', 'GC', 'GG', 'GT', 'AA', 'AG', 'AT', 'CA',
        'CG', 'CT', 'GA', 'GG', 'GT', 'AA', 'AC', 'AT', 'CA', 'CC', 'CT',
        'GA', 'GC', 'GT', 'AA', 'AC', 'AG', 'CA', 'CC', 'CG', 'GA', 'GC',
        'GG'
    ]
    order_m = ['{}_{}'.format(el[0], el[1]) for el in zip(types, subtypes)]

    types = np.asarray(types, dtype='object')
    subtypes = np.asarray(subtypes, dtype='object').T
    types_m = np.zeros((len(types), 1), dtype=object)
    subtypes_m = np.zeros((len(types), 1), dtype=object)
    for ix in range(len(types)):
        types_m[ix, 0] = types[ix]
        subtypes_m[ix, 0] = subtypes[ix]

    return types_m, subtypes_m, order_m


def indel_types():

    types = [
        'DEL_C', 'DEL_C', 'DEL_C', 'DEL_C', 'DEL_C', 'DEL_C',
        'DEL_T', 'DEL_T', 'DEL_T', 'DEL_T', 'DEL_T', 'DEL_T',
        'INS_C', 'INS_C', 'INS_C', 'INS_C', 'INS_C', 'INS_C',
        'INS_T', 'INS_T', 'INS_T', 'INS_T', 'INS_T', 'INS_T',
        'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats',
        'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats',
        'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats',
        'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats', 'DEL_repeats',
        'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats',
        'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats',
        'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats',
        'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats', 'INS_repeats',
        'DEL_MH', 'DEL_MH', 'DEL_MH', 'DEL_MH', 'DEL_MH', 'DEL_MH', 'DEL_MH', 'DEL_MH', 'DEL_MH', 'DEL_MH', 'DEL_MH'
    ]

    subtypes = [
        '1_1', '1_2', '1_3', '1_4', '1_5', '1_6+',
        '1_1', '1_2', '1_3', '1_4', '1_5', '1_6+',
        '1_0', '1_1', '1_2', '1_3', '1_4', '1_5+',
        '1_0', '1_1', '1_2', '1_3', '1_4', '1_5+',
        '2_1', '2_2', '2_3', '2_4', '2_5', '2_6+',
        '3_1', '3_2', '3_3', '3_4', '3_5', '3_6+',
        '4_1', '4_2', '4_3', '4_4', '4_5', '4_6+',
        '5+_1', '5+_2', '5+_3', '5+_4', '5+_5', '5+_6+',
        '2_0', '2_1', '2_2', '2_3', '2_4', '2_5+',
        '3_0', '3_1', '3_2', '3_3', '3_4', '3_5+',
        '4_0', '4_1', '4_2', '4_3', '4_4', '4_5+',
        '5+_0', '5+_1', '5+_2', '5+_3', '5+_4', '5+_5+',
        '2_1', '3_1', '3_2', '4_1', '4_2', '4_3', '5+_1', '5+_2', '5+_3', '5+_4', '5+_5+'
    ]

    order_m = ['{}_{}'.format(el[0], el[1]) for el in zip(types, subtypes)]
    types = np.asarray(types, dtype='object')
    subtypes = np.asarray(subtypes, dtype='object').T
    types_m = np.zeros((len(types), 1), dtype=object)
    subtypes_m = np.zeros((len(types), 1), dtype=object)
    for ix in range(len(types)):
        types_m[ix, 0] = types[ix]
        subtypes_m[ix, 0] = subtypes[ix]

    return types_m, subtypes_m, order_m

# pentamer + consequence
def snvs_signature_analyzer():

    nuc = ['A', 'C', 'G', 'T']
    midd = ['C', 'T']
    total = []

    # this is the reference
    for ref in midd:
        for mut in nuc:
            if mut != ref:
                for f in nuc:
                    for s in nuc:
                        for four in nuc:
                            for last in nuc:
                                a = '{}{}{}{}{}{}'.format(f, s, ref, four, last, mut)
                                total.append(a)
    return set(total)

def snvs_types():

    types = [
        'C>A', 'C>A', 'C>A', 'C>A', 'C>G', 'C>G', 'C>G', 'C>G', 'C>T',
        'C>T', 'C>T', 'C>T', 'T>A', 'T>A', 'T>A', 'T>A', 'T>C', 'T>C',
        'T>C', 'T>C', 'T>G', 'T>G', 'T>G', 'T>G', 'C>A', 'C>A', 'C>A',
        'C>A', 'C>G', 'C>G', 'C>G', 'C>G', 'C>T', 'C>T', 'C>T', 'C>T',
        'T>A', 'T>A', 'T>A', 'T>A', 'T>C', 'T>C', 'T>C', 'T>C', 'T>G',
        'T>G', 'T>G', 'T>G', 'C>A', 'C>A', 'C>A', 'C>A', 'C>G', 'C>G',
        'C>G', 'C>G', 'C>T', 'C>T', 'C>T', 'C>T', 'T>A', 'T>A', 'T>A',
        'T>A', 'T>C', 'T>C', 'T>C', 'T>C', 'T>G', 'T>G', 'T>G', 'T>G',
        'C>A', 'C>A', 'C>A', 'C>A', 'C>G', 'C>G', 'C>G', 'C>G', 'C>T',
        'C>T', 'C>T', 'C>T', 'T>A', 'T>A', 'T>A', 'T>A', 'T>C', 'T>C',
        'T>C', 'T>C', 'T>G', 'T>G', 'T>G', 'T>G'
    ]
    subtypes = [
        'ACA', 'ACC', 'ACG', 'ACT', 'ACA', 'ACC', 'ACG', 'ACT', 'ACA',
        'ACC', 'ACG', 'ACT', 'ATA', 'ATC', 'ATG', 'ATT', 'ATA', 'ATC',
        'ATG', 'ATT', 'ATA', 'ATC', 'ATG', 'ATT', 'CCA', 'CCC', 'CCG',
        'CCT', 'CCA', 'CCC', 'CCG', 'CCT', 'CCA', 'CCC', 'CCG', 'CCT',
        'CTA', 'CTC', 'CTG', 'CTT', 'CTA', 'CTC', 'CTG', 'CTT', 'CTA',
        'CTC', 'CTG', 'CTT', 'GCA', 'GCC', 'GCG', 'GCT', 'GCA', 'GCC',
        'GCG', 'GCT', 'GCA', 'GCC', 'GCG', 'GCT', 'GTA', 'GTC', 'GTG',
        'GTT', 'GTA', 'GTC', 'GTG', 'GTT', 'GTA', 'GTC', 'GTG', 'GTT',
        'TCA', 'TCC', 'TCG', 'TCT', 'TCA', 'TCC', 'TCG', 'TCT', 'TCA',
        'TCC', 'TCG', 'TCT', 'TTA', 'TTC', 'TTG', 'TTT', 'TTA', 'TTC',
        'TTG', 'TTT', 'TTA', 'TTC', 'TTG', 'TTT'
    ]

    order_m = ['{}_{}'.format(el[0], el[1]) for el in zip(types, subtypes)]
    types = np.asarray(types, dtype='object')
    subtypes = np.asarray(subtypes, dtype='object').T
    types_m = np.zeros((len(types), 1), dtype=object)
    subtypes_m = np.zeros((len(types), 1), dtype=object)
    for ix in range(len(types)):
        types_m[ix, 0] = types[ix]
        subtypes_m[ix, 0] = subtypes[ix]

    return types_m, subtypes_m, order_m


def get_types(type_mut):

    types_m = []
    subtypes_m = []
    order_m = []
    if type_mut == 'snv':
        types_m, subtypes_m, order_m = snvs_types()
    elif type_mut == 'indels':
        types_m, subtypes_m, order_m = indel_types()
    elif type_mut == 'dbs':
        types_m, subtypes_m, order_m = dbs_types()

    return types_m, subtypes_m, order_m


def order_muts(type_mut):

    order = []
    if type_mut == 'snv':
        first = ['A', 'C', 'G', 'T']
        pyr = ['C', 'T']
        for p in pyr:
            for mut in first:
                if mut != p:
                    for f in first:
                        for f2 in first:
                            comb = '{}[{}>{}]{}'.format(f, p, mut, f2)
                            order.append(comb)

    elif type_mut == 'mnvs':

        first = ['A', 'C', 'G', 'T']
        pyr = ['C', 'T']

        for p in pyr:
            for first_mutation in first:
                for second_mutation in first:
                    origin_triplet = '{}{}{}'.format(first_mutation, p, second_mutation)
                    for first_mutation in first:
                        for second_mutation in first:
                            for third_mutation in first:
                                last_triplet = '{}{}{}'.format(first_mutation, second_mutation, third_mutation)
                                order.append('{}_{}'.format(origin_triplet, last_triplet))
    else:
        types_m, subtypes_m, order = get_types(type_mut)

    return order


def order_to_plot_dbs():

    order = order_muts("dbs")
    listed = ['AC', 'AT', 'CC', 'CG', 'CT', 'GC', 'TA', 'TC', 'TG', 'TT']
    toplot = []
    for k in listed:
        for o in order:
            if o.split('_')[0] == k:
                toplot.append(o)

    return toplot


def order_to_plot_indel():

    order_plot = [
        'DEL_C_1_0', 'DEL_C_1_1', 'DEL_C_1_2', 'DEL_C_1_3', 'DEL_C_1_4', 'DEL_C_1_5+',
        'DEL_T_1_0', 'DEL_T_1_1', 'DEL_T_1_2', 'DEL_T_1_3', 'DEL_T_1_4', 'DEL_T_1_5+',
        'INS_C_1_0', 'INS_C_1_1', 'INS_C_1_2', 'INS_C_1_3', 'INS_C_1_4', 'INS_C_1_5+',
        'INS_T_1_0', 'INS_T_1_1', 'INS_T_1_2', 'INS_T_1_3', 'INS_T_1_4', 'INS_T_1_5+',
        'DEL_repeats_2_0', 'DEL_repeats_2_1', 'DEL_repeats_2_2', 'DEL_repeats_2_3', 'DEL_repeats_2_4',
        'DEL_repeats_2_5+', 'DEL_repeats_3_0', 'DEL_repeats_3_1', 'DEL_repeats_3_2', 'DEL_repeats_3_3',
        'DEL_repeats_3_4', 'DEL_repeats_3_5+', 'DEL_repeats_4_0', 'DEL_repeats_4_1', 'DEL_repeats_4_2',
        'DEL_repeats_4_3', 'DEL_repeats_4_4', 'DEL_repeats_4_5+', 'DEL_repeats_5+_0', 'DEL_repeats_5+_1',
        'DEL_repeats_5+_2', 'DEL_repeats_5+_3', 'DEL_repeats_5+_4', 'DEL_repeats_5+_5+',
        'INS_repeats_2_0', 'INS_repeats_2_1', 'INS_repeats_2_2', 'INS_repeats_2_3', 'INS_repeats_2_4',
        'INS_repeats_2_5+', 'INS_repeats_3_0', 'INS_repeats_3_1', 'INS_repeats_3_2', 'INS_repeats_3_3',
        'INS_repeats_3_4', 'INS_repeats_3_5+', 'INS_repeats_4_0', 'INS_repeats_4_1', 'INS_repeats_4_2',
        'INS_repeats_4_3', 'INS_repeats_4_4', 'INS_repeats_4_5+', 'INS_repeats_5+_0', 'INS_repeats_5+_1',
        'INS_repeats_5+_2', 'INS_repeats_5+_3', 'INS_repeats_5+_4', 'INS_repeats_5+_5+',
        'DEL_MH_2_1', 'DEL_MH_3_1', 'DEL_MH_3_2', 'DEL_MH_4_1', 'DEL_MH_4_2', 'DEL_MH_4_3',
        'DEL_MH_5+_1', 'DEL_MH_5+_2', 'DEL_MH_5+_3', 'DEL_MH_5+_4', 'DEL_MH_5+_5+',
    ]
    return order_plot


def order_to_plot_snvs():

    return order_muts('snv')


def order_sigprofiler_snvs():

    return sorted(order_muts('snv'))
