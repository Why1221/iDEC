#!/usr/bin/env bash

make 

echo "./idec -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -m 8 -f index"
./idec -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -m 8 -f index  
echo "./idec -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -f ./index"
./idec -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -f ./index
echo "./idec -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -k 5 -n 1192505 -q ./datasets/datasets/glove-hamming-128-test-compact.txt -r results/gnd-compact.txt"
./idec -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -k 5 -n 1192505 -q ./datasets/datasets/glove-hamming-128-test-compact.txt -r results/gnd-compact.txt  
echo "./idec -f index -k 5 -o results/idec.txt -r results/gnd-compact.txt -q ./datasets/datasets/glove-hamming-128-test-compact.txt -c 16.0 -es 0 -t 0.0005"
./idec -f index -k 5 -o results/idec.txt -r results/gnd-compact.txt -q ./datasets/datasets/glove-hamming-128-test-compact.txt -c 16.0 -es 0 -t 0.0005



