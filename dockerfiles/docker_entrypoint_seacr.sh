#!/bin/bash

# Make conda available
source /opt/conda/bin/activate

# Stack the envs so each toolset is on $PATH
conda activate cutrun
conda activate seacr
conda activate preseq
conda activate ucsc

# If the first argument ends in “.json”, hand it off to your pipeline script
if [[ "$1" =~ \.json$ ]]; then
  exec bash /usr/local/bin/cutrun.sh "$@"
else
  exec "$@"
fi
