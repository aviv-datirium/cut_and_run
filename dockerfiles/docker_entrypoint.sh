#!/bin/bash

# Automatically activate both environments
source /opt/conda/bin/activate
conda activate cutrun
conda activate macs2
conda activate preseq
conda activate ucsc

exec "$@"
