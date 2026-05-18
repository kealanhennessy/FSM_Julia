#!/usr/bin/env bash
# Regenerates every case's expected-*.tif oracle by invoking each
# case_*/run.sh. Requires FSM_EXE pointing at your own upstream C++
# fsm.exe build (this repo does not bundle the upstream C++). See the
# repository README, "Testing" -> Option B.
set -euo pipefail
cd "$(dirname "$0")"

if [[ -z "${FSM_EXE:-}" ]]; then
  echo "Error: FSM_EXE is not set." >&2
  echo "Point it at your upstream fsm.exe build, e.g.:" >&2
  echo "  FSM_EXE=/path/to/Barnes2020-FillSpillMerge/build/fsm.exe ./run_all.sh" >&2
  echo "See the repository README (Testing -> Option B) for build steps." >&2
  exit 1
fi

for case_dir in case_*/; do
  echo "=== $case_dir ==="
  (cd "$case_dir" && ./run.sh)
  echo
done

echo "All cases regenerated."
