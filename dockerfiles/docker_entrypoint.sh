#!/bin/bash

source /opt/conda/bin/activate
conda activate cutrun
conda activate macs2
conda activate preseq
conda activate ucsc

# If the first argument is a config file, run the pipeline
if [[ "$1" =~ \.json$ ]]; then
  exec bash /usr/local/bin/cutrun.sh "$@"
else
  exec "$@"
fi
