#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

for case_dir in case_*/; do
  echo "=== $case_dir ==="
  (cd "$case_dir" && ./run.sh)
  echo
done

echo "All cases regenerated."
