# import modules
import os
import sys

import click

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
scripts_path = os.path.abspath(os.path.join(__file__, '..'))
sys.path.insert(0, scripts_path)

from sign_utils.signature_processing import postprocessing_signatureanalyzer


# Global variables
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    '-i', '--input-signatures-path',
    required=True,
    help='Path of the signatures',
    type=click.Path()
)
@click.option(
    '-o', '--output-path',
    required=True,
    help='Output path where to save the results',
    type=click.Path()
)
def run(input_signatures_path, output_path):

    os.makedirs(output_path, exist_ok=True)

    exposures_types = [
        ('Pan_transcription.dbs', 'Pan_transcription.snvs', 'Pan_transcription.indels'),
        ('Pan_replication.dbs', 'Pan_replication.snvs', 'Pan_replication.indels'),
        ('Pan.dbs', 'Pan_full.snvs', 'Pan.indels')
    ]
    bias_types = ['Transcription', 'Replication', 'None']

    for exposures, bias in zip(exposures_types, bias_types):
        print("processing:", bias)
        for exposure in exposures:
            print('\tProcessing:', exposure.split('.')[-1])
            summaries = os.path.join(input_signatures_path, exposure, "{}.exposures.tsv".format(exposure))
            postprocessing_signatureanalyzer(summaries, output_path, bias=bias)


if __name__ == '__main__':
    run()
