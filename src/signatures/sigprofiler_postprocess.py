# Import modules
import os
import sys

import click

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
scripts_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.insert(0, scripts_path)

from sign_utils.signature_processing import create_summary_sigprofiler, \
    extract_best_signatures_sigprofiler, postprocessing_sigprofiler


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input-mat', required=True, help='Input matrix, without the suffixes', type=click.Path())
@click.option('-r', '--results-path', required=True, help='Path of the SigProfiler results', type=click.Path())
@click.option(
    '-o', '--output-path',
    required=True,
    help='Output directory where to save the results',
    type=click.Path()
)
def run(input_mat, results_path, output_path):

    os.makedirs(output_path, exist_ok=True)

    original_matrix_file_basename = os.path.basename(input_mat)
    original_matrix_file_snv = input_mat + '.snvs.dlm'
    original_matrix_file_dbs = input_mat + '.dbs.dlm'
    original_matrix_file_indels = input_mat + '.indels.dlm'

    for f in [original_matrix_file_snv, original_matrix_file_dbs, original_matrix_file_indels]:

        # create the summary of the extraction
        create_summary_sigprofiler(f, results_path)

        # plot how many active signatures we have to decide the number
        extract_best_signatures_sigprofiler(f, results_path, output_path, 0.8)

    # after looking at the plot created above, based on stability and
    # reconstruction error we decide the number of active signatures
    d_ttype_sigs = {
        original_matrix_file_basename + '.snvs': 25,
        original_matrix_file_basename + '.dbs': 5,
        original_matrix_file_basename + '.indels': 12,
    }

    postprocessing_sigprofiler(d_ttype_sigs, results_path, 0.8, output_path)


if __name__ == '__main__':
    run()
