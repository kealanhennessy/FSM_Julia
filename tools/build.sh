#!/usr/bin/env bash
# Build the dh_dump oracle utility.
#
# dh_dump.cpp uses the upstream C++ FSM headers (dephier + richdem) to dump
# the post-GetDepressionHierarchy state to a text file the Julia tests parse.
# This script compiles it standalone, without touching the upstream
# CMakeLists.
#
# Usage:
#   ./build.sh                              # uses ../../Barnes2020-FillSpillMerge as the C++ repo root
#   FSM_CPP=/path/to/repo ./build.sh        # override the C++ repo root
#
# Note: For the dumped depression vector to match Julia's
# get_depression_hierarchy bit-for-bit, the upstream
# `submodules/dephier/include/dephier/dephier.hpp` must have the
# deterministic outlet tie-break patch applied. See
# patches/dephier-deterministic-outlets.patch.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FSM_CPP="${FSM_CPP:-$SCRIPT_DIR/../../Barnes2020-FillSpillMerge}"

if [[ ! -d "$FSM_CPP" ]]; then
  echo "Error: FSM C++ repo not found at $FSM_CPP" >&2
  echo "Set FSM_CPP=/path/to/Barnes2020-FillSpillMerge and re-run." >&2
  exit 1
fi

DEPHIER_INC="$FSM_CPP/submodules/dephier/include"
RICHDEM_INC="$FSM_CPP/submodules/dephier/submodules/richdem/include"
FSM_INC="$FSM_CPP/include"

for d in "$DEPHIER_INC" "$RICHDEM_INC" "$FSM_INC"; do
  if [[ ! -d "$d" ]]; then
    echo "Error: include directory not found: $d" >&2
    echo "Run 'git submodule update --init --recursive' in $FSM_CPP." >&2
    exit 1
  fi
done

GDAL_CFLAGS="$(gdal-config --cflags)"
GDAL_LIBS="$(gdal-config --libs)"

# Homebrew's GDAL on this machine is built for x86_64; match it so the
# linker is happy regardless of the host architecture (Rosetta on
# Apple Silicon, native x86_64 otherwise). Override CXX/ARCH via env if
# you have a different setup.
CXX="${CXX:-/usr/local/opt/llvm/bin/clang++}"
ARCH="${ARCH:-x86_64}"

set -x
"$CXX" -arch "$ARCH" -std=c++17 -O2 \
  -DUSEGDAL \
  -I "$FSM_INC" \
  -I "$DEPHIER_INC" \
  -I "$RICHDEM_INC" \
  $GDAL_CFLAGS \
  "$SCRIPT_DIR/dh_dump.cpp" \
  $GDAL_LIBS \
  -o "$SCRIPT_DIR/dh_dump.exe"
