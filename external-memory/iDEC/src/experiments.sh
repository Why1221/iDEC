#!/usr/bin/env bash

echo "./srs -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train.txt -n 1192505 -m 8 -f index"
./srs -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train.txt -n 1192505 -m 8 -f index  
echo "./srs -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train.txt -n 1192505 -f ./index"
./srs -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train.txt -n 1192505 -f ./index 
echo "./srs -d 128 -ds ./datasets/datasets/glove-hamming-128-train.txt -k 5 -n 1192505 -q ./datasets/datasets/glove-hamming-128-test.txt -r results/gnd.txt"
./srs -d 128 -ds ./datasets/datasets/glove-hamming-128-train.txt -k 5 -n 1192505 -q ./datasets/datasets/glove-hamming-128-test.txt -r results/gnd.txt  
echo "./srs -f index -k 5 -o results/srs.txt -r results/gnd.txt -q ./datasets/datasets/glove-hamming-128-test.txt -c 16.0 -es 0 -t 0.0005"
./srs -f index -k 5 -o results/srs.txt -r results/gnd.txt -q ./datasets/datasets/glove-hamming-128-test.txt -c 16.0 -es 0 -t 0.0005