import pandas as pd
import subprocess
import sys
import os

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from config import Conf

def get_line(line):
    dic = {
        i.split()[0]: i.split()[1].replace('"', '') for i in line.split(";") if len(i.split()) > 1
    }
    return dic

# Cancer gene census
def get_CGC_genes():

    CGC = Conf['CGC_path']
    dfcgc = pd.read_csv(CGC, sep ='\t')
    append_syn = []
    for synonym in dfcgc['Synonyms']:
        if str(synonym)!= 'nan':
            for s in synonym.split(','):
                append_syn.append(s)
    cgc_list = set(dfcgc['Gene Symbol'].tolist()  + append_syn)
    return cgc_list

def exonic_regions():

    gencode_annotations = 'data/asymmetry_files/gencode.v29lift37.annotation.gtf.gz'

    gencode = pd.read_csv(
        gencode_annotations,
        compression='gzip',
        sep="\t",
        header=None,
        skiprows=5,
        usecols=[0, 2, 3, 4, 6, 8],
        names=['CHROMOSOME', 'TYPE', 'START', 'END', 'STRAND', 'INFOS']
    )

    gencode = gencode[~gencode['INFOS'].isnull()]

    # Remove the `chr` prefix from the CHROMOSOME column
    gencode['CHROMOSOME'] = gencode['CHROMOSOME'].map(lambda x: x[3:])

    # Parse the INFOS column
    gencode['INFOS_new'] = gencode['INFOS'].apply(get_line)

    gencode['GENE_ID'] = gencode['INFOS_new'].map(lambda x: x.get('gene_id', 'NaN.').split('.')[0])
    gencode['TRANSCRIPT_ID'] = gencode['INFOS_new'].map(lambda x: x.get('transcript_id', 'NaN.').split('.')[0])
    gencode['GENE_TYPE'] = gencode['INFOS_new'].map(lambda x: x.get('gene_type', 'NaN'))
    gencode['SYMBOL'] = gencode['INFOS_new'].map(lambda x: x.get('gene_name', 'NaN'))
    gencode['TRANSCRIPT_TYPE'] = gencode['INFOS_new'].map(lambda x: x.get('transcript_type', 'NaN'))

    gencode = gencode[gencode['START']!=gencode['END']]

    gencode = gencode[
        (gencode['GENE_TYPE'] == 'protein_coding') &
        (gencode['TRANSCRIPT_TYPE'] == 'protein_coding')
    ].copy()

    cds = gencode[(gencode['TYPE'] == 'CDS')].copy()

    cds = cds[['CHROMOSOME', 'START', 'END', 'STRAND', 'GENE_ID', 'GENE_ID', 'SYMBOL']].drop_duplicates().copy()
    cds['START'] = cds['START'].astype(int)
    cds['END'] = cds['END'].astype(int)
    cds.sort_values(by=['CHROMOSOME', 'START', 'END'], axis=0, ascending=True, inplace=True)
    cds.to_csv('data/asymmetry_files/hg19.cds.tsv.gz', sep="\t", compression='gzip', index=None, header=True)


    cmd = '''zcat data/asymmetry_files/hg19.cds.tsv.gz | tail -n +2 |awk '{OFS="\t"}{print $1, $2, $3, "1", "2",  $4, $4 "_" $5 "_" $6 "_" $7}' | bedtools merge -i stdin -s -c 7 -o distinct | sed 's/_/\t/g' | gzip > data/asymmetry_files/hg19.cds.merged.gz'''
    completed = subprocess.run(cmd, shell=True)

def get_cgc_exonic_regions():

    exonic_regions()
    cgc_list = get_CGC_genes()
    df = pd.read_csv('data/asymmetry_files/hg19.cds.merged.gz', sep = '\t',
                    names = ['CHROMOSOME', 'START', 'END', 'STRAND','GENE_ID', 'GENE_ID_2', 'SYMBOL'])

    genic_locations = df[df['SYMBOL'].isin(cgc_list)]

    genic_locations.to_csv('data/megabase_probability/cgc_exonic_regions.tsv',
                          sep ='\t', index = False, header = False)

if __name__ == '__main__':
    get_cgc_exonic_regions()
