#!/usr/bin/env bash
# Regenerates this case's expected-*.tif oracle from input.asc using the
# upstream C++ fsm.exe. This repo does not bundle the upstream C++, so
# there is no default path — point FSM_EXE at your own upstream build:
#
#   FSM_EXE=/path/to/Barnes2020-FillSpillMerge/build/fsm.exe ./run.sh
#
# See the repository README, "Testing" -> Option B, for how to produce
# that build.
set -euo pipefail
cd "$(dirname "$0")"

if [[ -z "${FSM_EXE:-}" ]]; then
  echo "Error: FSM_EXE is not set." >&2
  echo "Point it at your upstream fsm.exe build, e.g.:" >&2
  echo "  FSM_EXE=/path/to/Barnes2020-FillSpillMerge/build/fsm.exe ./run.sh" >&2
  echo "See the repository README (Testing -> Option B) for build steps." >&2
  exit 1
fi

gdal_translate -of GTiff -ot Float64 -q input.asc input.tif
"$FSM_EXE" input.tif expected -100 --swl 1.0
