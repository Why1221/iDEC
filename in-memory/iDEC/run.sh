#!/usr/bin/env bash

set -e

./idec -d 192 -n 54187 -ds ../EXPERIMENTS/audio/audio-hamming-192-train-compact.flat -if index -m 6  
./idec -d 192 -n 54187 -ds ../EXPERIMENTS/audio/audio-hamming-192-train-compact.flat -if index -qn 200 -qs ../EXPERIMENTS/audio/audio-hamming-192-test-compact.flat -rf idec-steal.txt -gt ../EXPERIMENTS/audio/gnd.txt -t 0.001
