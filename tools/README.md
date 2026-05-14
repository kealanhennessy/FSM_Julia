# tools/

Auxiliary scripts used during the Julia port for generating reference
data against the upstream C++ `Barnes2020-FillSpillMerge` algorithm.
None of this is needed at runtime — it only exists to regenerate the
oracle files under `test/test_cases/` when needed.

## Contents

- `dh_dump.cpp` — runs the C++ ocean-labeling + `GetDepressionHierarchy`
  on a GeoTIFF and writes the label grid, flowdirs grid, and full
  `Depression` vector to a plain-text file the Julia tests parse
  (`expected-dh.txt`). Built standalone via `build.sh` against the
  upstream C++ headers; the upstream `CMakeLists.txt` is not modified.
- `build.sh` — compiles `dh_dump.cpp` to `dh_dump.exe`. Set the env var
  `FSM_CPP` if your `Barnes2020-FillSpillMerge` checkout is not at
  `../../Barnes2020-FillSpillMerge` relative to this folder.
- `patches/` — local patches against the upstream C++ repo. These are
  required for `dh_dump.exe` (and `fsm.exe`) to produce output that
  matches the Julia port bit-for-bit, but they are **not** pushed
  upstream — the upstream repo belongs to a different author. Reapply
  them in any fresh checkout before rebuilding the oracle binaries:

  ```
  cd $FSM_CPP/submodules/dephier
  git apply $FSM_JULIA/tools/patches/dephier-deterministic-outlets.patch
  cd $FSM_CPP/submodules/dephier/submodules/richdem
  git apply $FSM_JULIA/tools/patches/richdem-gdal-cslconstlist.patch
  ```

## Regenerating oracle data

```
# 1. Apply patches (above) and build dh_dump
./build.sh

# 2. For each test case, dump the dephier output:
./dh_dump.exe ../test/test_cases/case_01_trough/input.tif 0.0 \
              ../test/test_cases/case_01_trough/expected-dh.txt

# The wtd oracle (expected-wtd.tif) is produced by fsm.exe from the
# upstream CMake build; rebuild it the usual way after applying the
# patches above.
```
