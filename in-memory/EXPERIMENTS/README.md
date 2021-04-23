# Experiments for In-Memory iDEC Algorithms

## Usage

+ Formating datasets
```bash
./prepare.py
```

+ Run experiments
```bash
Usage: run.sh [-ahr]

This script attempts to run ANN experiments
options:
 -a: (default) run (A)ll algorithms - good luck!
 -c: clean all!
 -h: print this (H)elp message
 -r <alg>: only run <alg>, where <alg>=LinearScan|iDEC
```
## Running Example
```bash
./run.sh -r LinearScan
./run.sh -r iDEC
```
Make sure to run **LinearScan** before running **iDEC** algorithms.
