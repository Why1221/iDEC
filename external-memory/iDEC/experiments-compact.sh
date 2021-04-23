#!/usr/bin/env bash

make idec-compact


echo "indexing ..."
echo "./idec-compact -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -m 6 -f index-for-compact"
./idec-compact -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -m 6 -f index-for-compact  

echo "create blocked binary data ..."
echo "./idec-compact -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -f ./index-for-compact"
./idec-compact -b 4096 -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -n 1192505 -f ./index-for-compact

echo "calculate ground truth ..."
echo "./idec-compact -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -k 1 -n 1192505 -q ./datasets/datasets/glove-hamming-128-test-compact.txt -r results/gnd-compact.txt"
./idec-compact -d 128 -ds ./datasets/datasets/glove-hamming-128-train-compact.bin -k 1 -n 1192505 -q ./datasets/datasets/glove-hamming-128-test-compact.txt -r results/gnd-compact.txt  

echo "query ..."
echo "./idec-compact -f index-for-compact -k 1 -o results/idec-compact.txt -r results/gnd-compact.txt -q ./datasets/datasets/glove-hamming-128-test-compact.txt -c 16.0 -es 0 -t 0.0005"
./idec-compact -f index-for-compact -k 1 -o results/idec-compact.txt -r results/gnd-compact.txt -q ./datasets/datasets/glove-hamming-128-test-compact.txt -c 16.0 -es 0 -t 0.0005



