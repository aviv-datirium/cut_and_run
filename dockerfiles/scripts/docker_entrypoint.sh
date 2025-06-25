#!/bin/bash
set -e

# Where your output lives; use /pipeline if that’s guaranteed to be mounted
TARGET_DIR=/mnt/data/home/aviv

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
