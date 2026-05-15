#!/usr/bin/env bash
# Build the C++ helper tools that feed FSM_Julia's test suite.
#
# Builds (always together — they're cheap):
#   dh_dump.exe              — Phase 2 dephier oracle dumper (text file
#                              for test/test_cases/case_*/expected-dh.txt)
#   gen_random_terrains.exe  — Phase 4 perlin terrain batch generator
#                              for test/test_cases/random/
#
# Neither tool touches the upstream CMakeLists; everything compiles
# standalone against a minimal vendored snapshot of the C++ FSM headers
# under ../vendor/Barnes2020-FillSpillMerge/. The two upstream patches
# (deterministic dephier tie-break + GDAL CSLConstList fix) are
# **pre-applied** in the vendored copy — no patch step needed before
# building. See the repository root README.md ("Testing" and "The
# vendored C++ snapshot" sections) for full details.
#
# Usage:
#   ./build.sh                              # uses ../vendor/Barnes2020-FillSpillMerge
#   FSM_CPP=/path/to/repo ./build.sh        # override (e.g. to use a full upstream clone)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FSM_CPP="${FSM_CPP:-$SCRIPT_DIR/../vendor/Barnes2020-FillSpillMerge}"

if [[ ! -d "$FSM_CPP" ]]; then
  echo "Error: FSM C++ repo not found at $FSM_CPP" >&2
  echo "Set FSM_CPP=/path/to/Barnes2020-FillSpillMerge and re-run." >&2
  exit 1
fi

DEPHIER_INC="$FSM_CPP/submodules/dephier/include"
RICHDEM_INC="$FSM_CPP/submodules/dephier/submodules/richdem/include"
FSM_INC="$FSM_CPP/include"
RICHDEM_SRC="$FSM_CPP/submodules/dephier/submodules/richdem/src/terrain_generation"

for d in "$DEPHIER_INC" "$RICHDEM_INC" "$FSM_INC"; do
  if [[ ! -d "$d" ]]; then
    echo "Error: include directory not found: $d" >&2
    if [[ "$FSM_CPP" == *"/vendor/"* ]]; then
      echo "The vendored snapshot appears incomplete. Re-clone FSM_Julia or repopulate vendor/." >&2
    else
      echo "Run 'git submodule update --init --recursive' in $FSM_CPP." >&2
    fi
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

CXX_COMMON=(
  "$CXX" -arch "$ARCH" -std=c++17 -O2
  -DUSEGDAL
  -I "$FSM_INC"
  -I "$DEPHIER_INC"
  -I "$RICHDEM_INC"
  $GDAL_CFLAGS
)

set -x
"${CXX_COMMON[@]}" \
  "$SCRIPT_DIR/dh_dump.cpp" \
  $GDAL_LIBS \
  -o "$SCRIPT_DIR/dh_dump.exe"

"${CXX_COMMON[@]}" \
  "$SCRIPT_DIR/gen_random_terrains.cpp" \
  "$RICHDEM_SRC/terrain_generation.cpp" \
  "$RICHDEM_SRC/PerlinNoise.cpp" \
  $GDAL_LIBS \
  -o "$SCRIPT_DIR/gen_random_terrains.exe"
