#!/bin/bash
# Make conda available
source /opt/conda/bin/activate

# Activate only the unified env
conda activate cutrun

# Dispatch to your pipeline or to whatever command the user passed
if [[ "$1" =~ \.json$ ]]; then
  exec bash /usr/local/bin/cutrun.sh "$@"
else
  exec "$@"
fi
