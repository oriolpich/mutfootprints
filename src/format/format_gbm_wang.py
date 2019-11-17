# Import modules
import pandas as pd


def format_data():

    f = 'data/GBM_Wang/ng.3590-S4.xlsx'
    somatic_muts = 'data/GBM_Wang/gbm_muts.tsv'

    print('Formating files...')
    df = pd.read_excel(f)
    df.to_csv(somatic_muts, sep='\t', header=True, index=False)
    df = pd.read_csv(somatic_muts, sep='\t')

    wantedchroms = [str(i) for i in range(1, 23)]
    wantedchroms.append('Y')
    wantedchroms.append('X')

    # select only variants with the PASS filter and in the chromosomes we are interested in
    df = df[df['Chr'].isin(wantedchroms)]

    # number of characters in alt
    df['len_alt'] = df['ref'].str.len()

    # number of characters in ref
    df['len_ref'] = df['alt'].str.len()

    # first classification between SNV and others
    df['TYPE'] = df.apply(
        lambda x: 'SNV' if (x['len_alt'] == 1) & (x['len_ref'] == 1) and
                           (x['alt'] != '-') and (x['ref'] != '-') else 'INDEL', axis=1
    )

    snvs = df[df['TYPE'] == 'SNV']

    snvs.sort_values(by='CaseID')[['CaseID', 'Chr', 'position', 'ref', 'alt']].to_csv(
        'data/GBM_Wang/gbm_muts.to_deconstruct.tsv', sep='\t', index=False, header=False
    )


if __name__ == '__main__':
    format_data()
