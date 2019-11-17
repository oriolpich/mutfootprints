import pandas as pd

def get_order_context():
    context = ["C.A.in.ACA","C.A.in.ACC","C.A.in.ACG","C.A.in.ACT","C.A.in.CCA",
               "C.A.in.CCC","C.A.in.CCG","C.A.in.CCT","C.A.in.GCA","C.A.in.GCC",
               "C.A.in.GCG","C.A.in.GCT","C.A.in.TCA","C.A.in.TCC","C.A.in.TCG",
               "C.A.in.TCT","C.G.in.ACA","C.G.in.ACC","C.G.in.ACG","C.G.in.ACT",
               "C.G.in.CCA","C.G.in.CCC","C.G.in.CCG","C.G.in.CCT","C.G.in.GCA",
               "C.G.in.GCC","C.G.in.GCG","C.G.in.GCT","C.G.in.TCA","C.G.in.TCC",
               "C.G.in.TCG","C.G.in.TCT","C.T.in.ACA","C.T.in.ACC","C.T.in.ACG",
               "C.T.in.ACT","C.T.in.CCA","C.T.in.CCC","C.T.in.CCG","C.T.in.CCT",
               "C.T.in.GCA","C.T.in.GCC","C.T.in.GCG","C.T.in.GCT","C.T.in.TCA",
               "C.T.in.TCC","C.T.in.TCG","C.T.in.TCT","T.A.in.ATA","T.A.in.ATC",
               "T.A.in.ATG","T.A.in.ATT","T.A.in.CTA","T.A.in.CTC","T.A.in.CTG",
               "T.A.in.CTT","T.A.in.GTA","T.A.in.GTC","T.A.in.GTG","T.A.in.GTT",
               "T.A.in.TTA","T.A.in.TTC","T.A.in.TTG","T.A.in.TTT","T.C.in.ATA",
               "T.C.in.ATC","T.C.in.ATG","T.C.in.ATT","T.C.in.CTA","T.C.in.CTC",
               "T.C.in.CTG","T.C.in.CTT","T.C.in.GTA","T.C.in.GTC","T.C.in.GTG",
               "T.C.in.GTT","T.C.in.TTA","T.C.in.TTC","T.C.in.TTG","T.C.in.TTT",
               "T.G.in.ATA","T.G.in.ATC","T.G.in.ATG","T.G.in.ATT","T.G.in.CTA",
               "T.G.in.CTC","T.G.in.CTG","T.G.in.CTT","T.G.in.GTA","T.G.in.GTC",
               "T.G.in.GTG","T.G.in.GTT","T.G.in.TTA","T.G.in.TTC","T.G.in.TTG","T.G.in.TTT"]

    order = [  ]
    first = ['A', 'C', 'G', 'T']
    pyr = ['C', 'T']
    for p in pyr:
        for mut in first:
            if mut != p:
                for f in first:
                    for f2 in first:
                        comb = '{}[{}>{}]{}'.format(f, p, mut, f2)
                        order.append(comb)
    return context, order

def prepare_HDP_input():
    context, order = get_order_context()

    path_matrix = 'data/hartwig/signatures/matrix/Colon-Rectum.snvs.dlm'
    outpath = 'data/hartwig/signatures/matrix/Colon-Rectum.HDP.txt'

    df = pd.read_csv(path_matrix, sep ='\t')
    df.index = sorted(order)
    df = df.T

    new_cols = ['{}.{}.in.{}{}{}'.format(col[2], col[4], col[0], col[2], col[-1]) for col in df.columns]
    df.columns = new_cols
    df = df[context]

    df.to_csv(outpath, sep='\t', index=True, header=True)

if __name__ == '__main__':
    prepare_HDP_input()

