#!/usr/bin/python3

import sys
import struct
import numpy as np
import h5py
import os.path

def get_enc_dim(h5filename):
    name, _ = os.path.splitext(h5filename)
    return int(int(name.split('-')[-1]) / 64)

def uint64_to_binary(num):
    return [c for c in f"{num:0>64b}"]

h5filename = 'datasets/glove-hamming-128.h5'
max_k = 100

enc_dim = get_enc_dim(h5filename)
with h5py.File(h5filename, 'r') as inf:
    with open('datasets/glove-hamming-128-train-compact.bin', 'wb') as oubcdf:
        with open('datasets/glove-hamming-128-train.txt', 'w') as oudf:
            with open('datasets/glove-hamming-128-test-compact.txt', 'w') as oucqf:
                with open('datasets/glove-hamming-128-test.txt', 'w') as ouqf:
                    counter = 1
                    train = inf['train'][:]
                    train.tofile(oubcdf)
                    n = (int)(train.shape[0] / enc_dim)
                    print(f"#dim: {enc_dim * 64}, #points: {n}")
                    for i in range(n):
                        bs = []
                        for j in range(enc_dim):
                            bs.extend(uint64_to_binary(train[i * enc_dim + j]))
                        sbs = " ".join(bs)
                        # it seems srs only accepts integer index
                        oudf.write(f'{i + 1} {sbs}\n')
                        
                        if counter % 10000 == 0:
                            sys.stdout.write('%d points processed...\n' % counter)
                        counter += 1
                    
                    query = inf['test']
                    qn = (int)(query.shape[0] / enc_dim)
                    for i in range(qn):
                        bs = []
                        for j in range(enc_dim):
                            bs.extend(uint64_to_binary(query[i * enc_dim + j]))
                        sbs = " ".join(bs)
                        ouqf.write(f'{i} {sbs}\n')
                        sbs = " ".join([str(q) for q in query[(i*enc_dim):((i+1)*enc_dim)]])
                        oucqf.write(f'{i} {sbs}\n')                     

