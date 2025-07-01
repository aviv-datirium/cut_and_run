#!/bin/bash
set -e

export PATH="/opt/conda/envs/cutrun/bin:$PATH"

# Dynamically determine the current working directory inside the container
TARGET_DIR=${TARGET_DIR:-$(pwd)}

uid=$(stat -c %u "$TARGET_DIR")
gid=$(stat -c %g "$TARGET_DIR")

# Create matching group & user only if needed
if ! getent group "$gid" >/dev/null; then
  groupadd -g "$gid" hostgrp
fi
if ! id -u "$uid" >/dev/null 2>&1; then
  useradd -m -u "$uid" -g "$gid" hostusr
fi

exec gosu "$uid:$gid" micromamba run -n cutrun /usr/local/bin/cutrun.sh "$@"
