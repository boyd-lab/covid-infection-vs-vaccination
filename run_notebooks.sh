#!/bin/bash
# Run a set of notebooks
# Optionally set kernel with --set-kernel "kernel-name"

set -euo pipefail

for notebook in "$@"
do
    echo "Running $notebook";
    jupytext --sync --pipe black --execute "$notebook";
    echo;
done
