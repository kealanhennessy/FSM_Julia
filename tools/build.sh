#!/usr/bin/env bash
# Build the C++ helper tools that feed FSM_Julia's test suite.
#
# Builds (always together — they're cheap):
#   dh_dump.exe              — dephier oracle dumper (text file for
#                              test/test_cases/case_*/expected-dh.txt)
#   pf_dump.exe              — Zhou2016 Priority-Flood oracle dumper
#                              (filled DEM .tif for test/test_cases/random_pf/)
#   gen_random_terrains.exe  — perlin terrain batch generator for
#                              test/test_cases/random/
#
# Neither tool touches the upstream CMakeLists; everything compiles
# standalone against the headers in your own clone of the upstream C++
# Barnes2020-FillSpillMerge. There is no bundled copy — you point
# FSM_CPP at a clone that you have checked out at the pinned commit and
# patched. See the repository root README.md, "Testing" -> Option B,
# for the exact clone/checkout/patch recipe (including the two required
# patches in tools/patches/ and the upstream commit SHAs).
#
# Usage:
#   FSM_CPP=/path/to/Barnes2020-FillSpillMerge ./build.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "${FSM_CPP:-}" ]]; then
  echo "Error: FSM_CPP is not set." >&2
  echo "Point it at your patched upstream Barnes2020-FillSpillMerge clone:" >&2
  echo "  FSM_CPP=/path/to/Barnes2020-FillSpillMerge ./build.sh" >&2
  echo "See the repository README (Testing -> Option B) for the clone +" >&2
  echo "checkout + patch recipe." >&2
  exit 1
fi

if [[ ! -d "$FSM_CPP" ]]; then
  echo "Error: FSM_CPP=$FSM_CPP is not a directory." >&2
  exit 1
fi

DEPHIER_INC="$FSM_CPP/submodules/dephier/include"
RICHDEM_INC="$FSM_CPP/submodules/dephier/submodules/richdem/include"
FSM_INC="$FSM_CPP/include"
RICHDEM_SRC="$FSM_CPP/submodules/dephier/submodules/richdem/src/terrain_generation"

for d in "$DEPHIER_INC" "$RICHDEM_INC" "$FSM_INC"; do
  if [[ ! -d "$d" ]]; then
    echo "Error: include directory not found: $d" >&2
    echo "Run 'git submodule update --init --recursive' in $FSM_CPP" >&2
    echo "(and check it out at the pinned commit — see README Option B)." >&2
    exit 1
  fi
done

GDAL_CFLAGS="$(gdal-config --cflags)"
GDAL_LIBS="$(gdal-config --libs)"

# Compiler: default to the system C++ compiler on PATH (Apple clang on
# macOS, g++/clang++ on Linux). Override CXX for a specific toolchain,
# e.g. Homebrew LLVM:  CXX=/usr/local/opt/llvm/bin/clang++ ./build.sh
CXX="${CXX:-c++}"

# Target architecture: by default, build for the host's native arch
# (don't pass -arch at all). Set ARCH only if your GDAL was built for a
# different architecture than your host. The classic case: Apple
# Silicon with an x86_64 Homebrew GDAL (Homebrew installed under
# /usr/local via Rosetta) — there you need:  ARCH=x86_64 ./build.sh
# (the compiler cross-targets; the binary then runs under Rosetta).
# The Julia tests don't care about the helper binaries' architecture —
# they only consume the data the binaries produce.
CXX_COMMON=( "$CXX" )
if [[ -n "${ARCH:-}" ]]; then
  CXX_COMMON+=( -arch "$ARCH" )
fi
CXX_COMMON+=(
  -std=c++17 -O2
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
  "$SCRIPT_DIR/pf_dump.cpp" \
  $GDAL_LIBS \
  -o "$SCRIPT_DIR/pf_dump.exe"

"${CXX_COMMON[@]}" \
  "$SCRIPT_DIR/gen_random_terrains.cpp" \
  "$RICHDEM_SRC/terrain_generation.cpp" \
  "$RICHDEM_SRC/PerlinNoise.cpp" \
  $GDAL_LIBS \
  -o "$SCRIPT_DIR/gen_random_terrains.exe"
