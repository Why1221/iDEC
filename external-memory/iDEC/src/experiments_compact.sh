#!/usr/bin/env bash

echo "./srs-compact -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -m 8 -f index-for-compact"
./srs-compact -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -m 8 -f index-for-compact  
echo "./srs-compact -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -f ./index-for-compact"
./srs-compact -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -f ./index-for-compact 
echo "./srs-compact -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -k 5 -n 1192505 -q ./datasets/datasets/glove-hamming-128-test-compact.txt -r results/gnd-compact.txt"
./srs-compact -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -k 5 -n 1192505 -q ./datasets/datasets/glove-hamming-128-test-compact.txt -r results/gnd-compact.txt  
echo "./srs-compact -f index-for-compact -k 5 -o results/srs-compact.txt -r results/gnd-compact.txt -q ./datasets/datasets/glove-hamming-128-test-compact.txt -c 16.0 -es 0 -t 0.0005"
./srs-compact -f index-for-compact -k 5 -o results/srs-compact.txt -r results/gnd-compact.txt -q ./datasets/datasets/glove-hamming-128-test-compact.txt -c 16.0 -es 0 -t 0.0005



