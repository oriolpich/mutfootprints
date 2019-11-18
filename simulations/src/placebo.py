import os
import json
import click
import gzip
import pickle
import numpy as np
import tqdm

from cotreatment import attribution, attribution_many_to_one


@click.group()
def cli():
    pass


@cli.command()
@click.option('--signatures_path', type=click.Path())
@click.option('--deconstruction_folder', type=click.Path())
@click.option('--matrix_treatment_folder', type=click.Path())
@click.option('--output_folder', type=click.Path())
@click.option('--symmetric', is_flag=True)
@click.option('--one_to_one', is_flag=True)
def run_placebo(signatures_path, deconstruction_folder, matrix_treatment_folder, output_folder, symmetric, one_to_one):

    with gzip.open(signatures_path, 'rb') as f:
        signatures = pickle.load(f)

    for overlap in tqdm.tqdm([25, 50, 75]):
        for sig in tqdm.tqdm(['SBS9', 'SBS31']):
            global_pool = []
            for n_treated in [10, 25, 50, 100, 150]:
                replicate_pool = []
                for replicate in range(1, 11):
                    try:
                        if one_to_one:
                            res = attribution(sig, n_treated, replicate,
                                              deconstruction_folder, matrix_treatment_folder,
                                              overlap, signatures, symmetric)
                        else:
                            res = attribution_many_to_one(sig, n_treated, replicate,
                                                          deconstruction_folder, matrix_treatment_folder,
                                                          overlap, signatures, symmetric)
                        explained_by_sig = np.array(res[sig]) + 1
                        explained_by_placebo = np.array(res['placebo']) + 1
                        fold_change = list(np.log10(explained_by_sig / explained_by_placebo))
                        replicate_pool += fold_change
                    except FileNotFoundError:
                        continue
                global_pool.append(replicate_pool)

            symmetric_label = 'symmetric' if symmetric else 'asymmetric'
            one_to_one_label = 'one_to_one' if one_to_one else 'many_to_one'
            fn = os.path.join(output_folder,
                         f'overlap.{overlap}.{sig}.{symmetric_label}.{one_to_one_label}.cotreatment_fold_change.json')

            with open(fn, 'wt') as f:
                json.dump(global_pool, f)


if __name__ == '__main__':

    cli()
