#!/usr/bin/env bash
set -e

VENV_DIR=".venv"

# Create venv if it doesn't exist
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment..."
    python -m venv $VENV_DIR
    $VENV_DIR/bin/pip install --upgrade pip
    $VENV_DIR/bin/pip install -e .
fi

# Always run inside the venv
$VENV_DIR/bin/python -m sgevalviz.cli "$@"