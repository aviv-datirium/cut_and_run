#!/bin/bash
set -e

export PATH="/opt/conda/envs/cutrun/bin:$PATH"

/usr/local/bin/cutrun.sh "$@"
