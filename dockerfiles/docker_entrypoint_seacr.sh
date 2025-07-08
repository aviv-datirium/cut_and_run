#!/bin/bash
set -e

# load conda’s shell function
source /opt/conda/etc/profile.d/conda.sh

# always activate the unified 'cutrun' environment
conda activate cutrun

# if the first argument is a .json config, run the pipeline script
if [[ "$1" =~ \.json$ ]]; then
  exec bash /usr/local/bin/cutrun.sh "$1"
else
  # otherwise, run whatever command was passed
  exec "$@"
fi
