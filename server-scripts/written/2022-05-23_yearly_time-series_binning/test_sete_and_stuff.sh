#!/bin/bash

set -ex

echo 0

source /home/rrohwer/miniconda3/bin/activate metabat2_robin

echo 1
echo 2
echo 3

metabat2 --help

echo 4

bbmap.sh -h

echo 5
