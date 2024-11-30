#!/bin/bash
set -e

# Download ViennaRNA tarball
wget https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.7.0/ViennaRNA-2.7.0.tar.gz -O /tmp/ViennaRNA-2.7.0.tar.gz

# Extract the tarball
tar -zxvf /tmp/ViennaRNA-2.7.0.tar.gz -C /tmp

# Navigate to the directory and install
cd /tmp/ViennaRNA-2.7.0
./configure
make
make install

# Cleanup
rm -rf /tmp/ViennaRNA-2.7.0 /tmp/ViennaRNA-2.7.0.tar.gz
