#!/usr/bin/env bash
# Make a self-contained folder for Shiny Server: app.py, requirements.txt, README.md, LICENSE and a copy of the package.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:?usage: make_bundle.sh OUTPUT_FOLDER}"
mkdir -p "$out"
cp "$here/app.py" "$here/requirements.txt" "$here/README.md" "$out/"
cp "$here/../../LICENSE" "$out/LICENSE"
rsync -a --delete --exclude '__pycache__' "$here/../barriernet/" "$out/barriernet/"
echo "bundle written to $out"
