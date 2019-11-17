# Import modules
import os
import click


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def write_chromosome(output_folder, chromosome, sequence):
    print('Chromosome: {}'.format(chromosome))
    with open(os.path.join(output_folder, '{}.txt'.format(chromosome)), 'w') as fo:
        fo.write(sequence)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input-fasta', required=True, help='Input fasta file', type=click.Path())
@click.option('-o', '--output-folder', required=True, help='Output directory', type=click.Path())
def run(input_fasta, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    with open(input_fasta, 'r') as fd:
        sequence = ''
        chromosome = ''

        for line in fd:
            if line.startswith(">"):
                if sequence != '':
                    write_chromosome(output_folder, chromosome, sequence)
                chromosome = line[1:].split()[0]
                sequence = ''
                continue
            sequence += line.strip()
        if sequence != '':
            write_chromosome(output_folder, chromosome, sequence)


if __name__ == '__main__':
        run()
