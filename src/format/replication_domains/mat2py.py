# Import modules
import click
import pandas as pd
import scipy.io


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input-mat', required=True, help='Input matrix', type=click.Path())
@click.option('-o', '--output-bed', required=True, help='Output file in bed format', type=click.Path())
def run(input_mat, output_bed):
    mat = scipy.io.loadmat(input_mat)
    data = mat['W'][0][0]
    names = mat['W'][0].dtype.names
    columns = ['chr', 'st', 'en', 'is_left', 'is_right', 'txplus', 'txminus', 'rt']
    df = pd.concat(
        [pd.DataFrame(v, columns=[k]) for k, v in zip(names, data) if k in columns],
        axis=1
    )
    df['st'] -= 1
    df['chr'] = df['chr'].map(lambda x: 'chr' + str(x))
    df.sort_values(by=['chr', 'st', 'en'], inplace=True)
    df.fillna(value=0.0, inplace=True)
    df[columns].to_csv(
        output_bed, sep='\t', index=False, header=True, compression='gzip'
    )


if __name__ == '__main__':
    run()
