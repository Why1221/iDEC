#!/usr/bin/env bash

set -e 


cd /tmp
git clone https://github.com/google/googletest.git
cd googletest
mkdir mybuild       # Create a directory to hold the build output.
cd mybuild
cmake ..  # Generate native build scripts.
make && sudo make install