#!/bin/bash
# Usage: ./run_notebook_to_new_file.sh [INPUT_NOTEBOOK_FILE] [GENERATED_OUTPUT_FILE]
# Run a notebook and save out as a new notebook. Useful when executing a template
# Optionally set kernel with --set-kernel "kernel-name"
jupytext --sync --pipe black --execute --output "$2" "$1";
