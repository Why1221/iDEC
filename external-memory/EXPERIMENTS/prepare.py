#!/usr/bin/env python3

from enum import Enum
import sys
import struct
import numpy as np
import h5py
import os.path
import click


class Scale(Enum):
    SMALL = 0
    MEDIUM = 1
    LARGE = 2


DATASETS = {
    Scale.SMALL: ['audio', 'mnist', 'enron'],
    Scale.MEDIUM: ['glove', 'GIST1M', 'SIFT1M'],
    Scale.LARGE: ['GIST80M', "SIFT1B"]
}
DATASET_INFO = {}
CUR_PATH = os.path.dirname(__file__)


def get_enc_dim(h5filename):
    name, _ = os.path.splitext(h5filename)
    return int(int(name.split('-')[-1]) / 64)


def uint64_to_binary(num):
    return [c for c in f"{num:0>64b}"]


def get_dsname(h5filename):
    dsname, _, _ = h5filename.partition('-')
    return dsname


def convert(h5filename, odir):
    assert os.path.exists(odir)
    print(f'Converting {h5filename} ...')
    _, name = os.path.split(h5filename)
    filename_prefix, _ = os.path.splitext(name)
    # dsname = get_dsname(filename_prefix)
    enc_dim = get_enc_dim(h5filename)
    qn = -1
    with h5py.File(h5filename, 'r') as inf:
        train_fn = f"{odir}/{filename_prefix}-train-compact.bin"
        test_tfn = f"{odir}/{filename_prefix}-test-compact.txt"
        test_bfn = f"{odir}/{filename_prefix}-test-compact.bin"

        with open(train_fn, 'wb') as oubcdf:
            with open(test_bfn, 'wb') as oubqf:
                with open(test_tfn, 'w') as ouqf:
                    train = inf['train'][:]
                    train.tofile(oubcdf)
                    n = (int)(train.shape[0] / enc_dim)
                    print(f"#dim: {enc_dim * 64}, #points: {n}")

                    query = inf['test'][:]
                    query.tofile(oubqf)
                    qn = (int)(query.shape[0])
                    qn = (int)(query.shape[0] / enc_dim)
                    for i in range(qn):
                        sbs = " ".join(
                            [str(q) for q in query[(i*enc_dim):((i+1)*enc_dim)]])
                        ouqf.write(f'{i} {sbs}\n')
        return {
            'train-b': os.path.split(train_fn)[-1],
            'test-b': os.path.split(test_bfn)[-1],
            'test-t': os.path.split(test_tfn)[-1],
            'n': int(n),
            'dimension': int(enc_dim * 64),
            'qn': int(qn),
            'dsh5': name
        }


def prepare_datasets():
    ddir = os.path.join(CUR_PATH, 'datasets')
    for f in os.listdir(ddir):
        if f.endswith('.h5'):
            dsname = get_dsname(f)
            if os.path.exists(dsname):
                print("directory already exists")
                exit(1)
            os.makedirs(dsname)
            DATASET_INFO[dsname] = convert(os.path.join(ddir,f), dsname)
    with open('datasets_info.json', 'w') as jf:
        import json 
        json.dump(DATASET_INFO, jf, indent=4)
    with open('dataset_info.txt', 'w') as tf:
        # for ds, info in DATASET_INFO.items():
        for ds in sorted(DATASET_INFO.keys(), key = lambda k: (DATASET_INFO[k]['n'], DATASET_INFO[k]['dimension'])):
            info = DATASET_INFO[ds]
            trainb = info['train-b']
            testb = info['test-b']
            testt = info['test-t']
            n = info['n']
            dim = info['dimension']
            qn = info['qn']
            dsh5 = info['dsh5']
            tf.write(f"{ds}%{trainb}%{testb}%{testt}%{n}%{dim}%{qn}%{dsh5}\n")


def prepare_large_datasets():
    ddir = os.path.join(CUR_PATH, 'large-datasets')
    for f in os.listdir(ddir):
        if f.endswith('.h5'):
            dsname = get_dsname(f)
            if os.path.exists(dsname):
                print("directory already exists")
                exit(1)
            os.makedirs(dsname)
            DATASET_INFO[dsname] = convert(os.path.join(ddir,f), dsname)
    with open('datasets_info.json', 'w') as jf:
        import json 
        json.dump(DATASET_INFO, jf, indent=4)
    with open('dataset_info.txt', 'w') as tf:
        # for ds, info in DATASET_INFO.items():
        for ds in sorted(DATASET_INFO.keys(), key = lambda k: (DATASET_INFO[k]['n'], DATASET_INFO[k]['dimension'])):
            info = DATASET_INFO[ds]
            trainb = info['train-b']
            testb = info['test-b']
            testt = info['test-t']
            n = info['n']
            dim = info['dimension']
            qn = info['qn']
            dsh5 = info['dsh5']
            tf.write(f"{ds}%{trainb}%{testb}%{testt}%{n}%{dim}%{qn}%{dsh5}\n")

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--large/--no-large', default=False)
def prepare(large):
    """Prepare datasets."""
    prepare_datasets()
    if large:
        prepare_large_datasets()

if __name__ == '__main__':
    prepare()

