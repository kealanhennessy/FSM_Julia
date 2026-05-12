#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

FSM_EXE="${FSM_EXE:-/Users/kealan/Desktop/Barnes2020-FillSpillMerge/build/fsm.exe}"

gdal_translate -of GTiff -ot Float64 -q input.asc input.tif
"$FSM_EXE" input.tif expected 0 --swl 1.0
