# FSM_Julia

A Julia port of the C++ Fill-Spill-Merge (FSM) algorithm of
Barnes ([2020](https://github.com/r-barnes/Barnes2020-FillSpillMerge)).
The port is line-by-line faithful to the upstream so that outputs match
the reference `fsm.exe` bit-for-bit on the committed test cases.

This branch is **self-contained**: a minimal copy of the upstream C++
sources is vendored at `vendor/Barnes2020-FillSpillMerge/`, so cloning
the repo gives you everything needed to run the full test suite — no
external repositories to fetch, no patches to apply.

## Quick start

```sh
git clone <this repo>
cd FSM_Julia
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Expected: `Pkg.test()` reports `1665 / 1665 pass` in roughly 6 seconds.

That's the entire test flow. Nothing else needs to be built — every
oracle file the tests compare against lives in
`test/test_cases/`. The vendored C++ and the helpers in `tools/` are
only required if you want to *regenerate* those oracles.

## Step-by-step

### 1. Prerequisites

- **Julia** ≥ 1.11. The dev machine runs 1.12.6.
- A working network connection for `Pkg.instantiate()` to fetch the
  Julia dependencies (`ArchGDAL`, `ArgParse`, `DataStructures`, `Test`).
  Total install footprint: ~140 packages, a few hundred MB on first
  run. Subsequent runs use the local Julia depot cache.

You do **not** need a C++ toolchain, GDAL, or anything else to run the
tests. All test data is committed under `test/test_cases/` (~3.5 MB
total: 5 hand-crafted oracle cases for the depression-hierarchy and
WTD outputs, plus 70 perlin-generated terrains for property tests).

### 2. Clone and enter the repo

```sh
git clone <repo-url> FSM_Julia
cd FSM_Julia
```

### 3. Instantiate the Julia environment

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This reads `Project.toml` / `Manifest.toml` and downloads the exact
dependency versions used by the dev machine. Takes a minute or two on
first run.

### 4. Run the full test suite

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected output (final two lines):

```
Test Summary:  | Pass  Total  Time
FillSpillMerge | 1665   1665  ~6s
```

The breakdown:

| Group | Tests | What it checks |
|---|---|---|
| Phase 1 unit tests | ~30 | constants, fp-compare helpers, BucketFillFromEdges |
| Phase 2 helpers | ~30 | DisjointDenseIntSet, LIFOMinPriorityQueue, Depression / Outlet |
| Phase 2 dephier oracles | ~250 | bit-exact match against `expected-dh.txt` (5 cases) |
| Phase 3 WTD oracles | 15 | bit-exact match against `expected-wtd.tif` (5 cases, tolerance 1e-9 / 1e-6) |
| Phase 4 unit ports | ~32 | direct ports of every `TEST_CASE` in `unittests/fsm_tests.cpp` |
| Phase 4 random property tests | ~1130 | 5 properties iterated over 70 perlin terrains, including a heavy-flooding-vs-Priority-Flood cross-check |
| Phase 4 edge cases | ~70 | 1×1 grids, all-ocean, monotonic slope, deep pit, Float32 elevations |

### 5. (Optional) Regenerate oracle data

Only needed if you are refreshing the test fixtures (e.g. after a
behavioural change to the algorithm). Requires the vendored C++
snapshot to be rebuilt locally, which means you need a C++17 toolchain
and GDAL. See `tools/README.md` for full instructions; the short
version on a machine like the dev machine:

```sh
brew install llvm gdal       # if you don't already have them
tools/build.sh               # builds tools/dh_dump.exe + tools/gen_random_terrains.exe

# Regenerate the 5 dephier oracles
for d in test/test_cases/case_*; do
  tools/dh_dump.exe "$d/input.tif" 0.0 "$d/expected-dh.txt"
done

# Regenerate the 70-terrain random batch
tools/gen_random_terrains.exe test/test_cases/random
```

Regenerating `expected-wtd.tif` is **not** supported against the
minimal vendored copy — that file is produced by upstream's full
`fsm.exe`, which the vendor doesn't include. To refresh those, clone
the full upstream Barnes2020-FillSpillMerge, apply the two patches in
`tools/patches/`, and build via CMake.

## Reference system

These specs are the only configuration the port has been exercised on
end-to-end. Anything reasonably close should work, but the exact
numbers are recorded here so reproducibility issues can be triaged
against a known-good baseline.

| Component | Version |
|---|---|
| Hardware | Apple M1 Pro |
| RAM | 16 GB |
| OS | macOS Tahoe 26.3.1 |
| Julia | 1.12.6 |
| Homebrew LLVM clang | 22.1.5 (`/usr/local/opt/llvm/bin/clang++`, x86_64 prefix) |
| GDAL | 3.13.0 (Homebrew, x86_64) |

Note on the x86_64 prefix: `tools/build.sh` defaults to `CXX=/usr/local/opt/llvm/bin/clang++` and `-arch x86_64` because the dev machine's Homebrew is installed under `/usr/local` (the Intel prefix, running through Rosetta on Apple Silicon) and its GDAL is an x86_64 build. Override either via env var:

```sh
CXX=clang++ ARCH=arm64 tools/build.sh
```

if you have a native ARM toolchain + GDAL. The Julia tests don't care
about the helper binaries' architecture (they only consume the data
the binaries produce).

## Repository layout

```
src/                          Julia port (this is the algorithm)
  FillSpillMerge.jl           Top-level module; just includes the rest
  constants.jl                D8 offsets, OCEAN / NO_DEP sentinels
  fp_compare.jl               fp_eq / fp_le / fp_ge (1e-6 tolerance)
  bucket_fill.jl              BucketFillFromEdges (Phase 1)
  types.jl                    Depression / Outlet structs (Phase 2)
  disjoint_set.jl             DisjointDenseIntSet (Phase 2)
  priority_queue.jl           LIFOMinPriorityQueue (Phase 2)
  dephier.jl                  GetDepressionHierarchy (Phase 2)
  fill_spill_merge.jl         FillSpillMerge driver + 9 helpers (Phase 3)
  priority_flood.jl           Zhou2016 Priority-Flood for cross-check (Phase 4)

test/
  runtests.jl                 Test runner entry point
  dephier_oracle.jl           Phase 2 oracle helpers (expected-dh.txt parser)
  fsm_unit_tests.jl           Phase 4 deterministic C++ port (fsm_tests.cpp)
  fsm_random_tests.jl         Phase 4 random property tests
  fsm_edge_tests.jl           Phase 4 edge cases + Float32
  test_cases/
    case_01_trough/ ... case_05_perlin_100x100/   Hand-crafted oracle cases
    random/                                         70 perlin terrains for property tests

tools/                        C++ helpers + patches (only for regenerating oracles)
  dh_dump.cpp                 Dumps GetDepressionHierarchy to expected-dh.txt
  gen_random_terrains.cpp     Generates the 70-terrain batch
  build.sh                    Builds both helpers against vendor/
  patches/                    For-reference diffs vs. upstream (pre-applied in vendor)
  README.md                   Details on regenerating oracle data

vendor/Barnes2020-FillSpillMerge/   Minimal vendored upstream C++ snapshot
  include/                          FSM C++ header
  submodules/dephier/include/       dephier headers (deterministic-outlets patch pre-applied)
  submodules/dephier/submodules/
    richdem/include/                richdem headers (GDAL-CSLConstList patch pre-applied)
    richdem/src/terrain_generation/ perlin source (linked by gen_random_terrains)
  LICENSE                           MIT (FSM)
  submodules/dephier/LICENSE        MIT (dephier)
  submodules/dephier/submodules/richdem/LICENSE.txt   GPL v3 (richdem)
  README.md                         Vendor scope, upstream commits, patch details
```

## Project phases

The port was developed across five distinct phases. Each phase is a
discrete chunk of functionality with its own oracle, so regressions
trace back to the phase that introduced them.

| Phase | Adds | Oracle |
|---|---|---|
| 0 | C++ build, 5 hand-crafted test cases | n/a (setup) |
| 1 | Constants, fp-compare, BucketFillFromEdges | unit tests |
| 2 | DisjointDenseIntSet, priority queue, dephier port | `expected-dh.txt` (bit-exact, via patched upstream) |
| 3 | `fill_spill_merge.jl` (~430 lines, all 10 functions) | `expected-wtd.tif` (1e-9, 1e-6 for case_05) |
| 4 | Zhou2016 Priority-Flood port, full C++ unit-test port (deterministic + 5 randomized + edge cases + Float32) | direct ports of `unittests/fsm_tests.cpp` + property tests |

## Licenses

The vendored C++ subtree at `vendor/Barnes2020-FillSpillMerge/`
preserves the original upstream licenses:

- FSM (top-level): MIT, Copyright (c) 2020 Richard Barnes
- dephier submodule: MIT, Copyright (c) 2020 Richard Barnes
- richdem submodule: **GPL v3** — note that the vendored richdem
  subdirectory carries GPL3 obligations (preserve `LICENSE.txt`, offer
  source on redistribution).

The Julia code under `src/` does not link against richdem at runtime;
it only reads files produced by separately-invoked C++ binaries
(`tools/*.exe`), so the GPL's "mere aggregation" reading applies. The
Julia port's own licensing is left to a top-level `LICENSE` file (not
yet set on this branch).
