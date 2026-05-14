# tools/

Auxiliary tooling used to (re)generate the reference data the Julia
tests run against. None of this runs as part of `Pkg.test()` — the
committed oracle files under `test/test_cases/` are what the tests
actually consume. These tools only matter when you want to refresh
those oracles.

Everything builds against the **minimal vendored snapshot** of the
upstream `Barnes2020-FillSpillMerge` C++ project at
`../vendor/Barnes2020-FillSpillMerge/`. The two required upstream
patches are **pre-applied** in the vendored copy, so a cold clone of
FSM_Julia can run `./build.sh` immediately. See
`../vendor/Barnes2020-FillSpillMerge/README.md` for the vendor scope,
source commits, and patch details.

## Contents

- `dh_dump.cpp` — runs the C++ ocean-labeling + `GetDepressionHierarchy`
  on a GeoTIFF and writes the label grid, flowdirs grid, and full
  `Depression` vector to a plain-text file the Phase 2 tests parse
  (`expected-dh.txt` under each `test/test_cases/case_*/`).
- `gen_random_terrains.cpp` — uses richdem's `perlin` to generate a
  fixed batch of small + large random terrains under
  `test/test_cases/random/`. Used by the Phase 4 property tests
  (mass conservation, repeated FSM, etc.). Deterministic — bumping
  `MASTER_SEED` in the source produces a fresh batch.
- `build.sh` — compiles both `.cpp` files to `.exe` binaries. Defaults
  to the vendored C++ snapshot at `../vendor/Barnes2020-FillSpillMerge`.
  Override the C++ root by setting the `FSM_CPP` env var (e.g. point
  it at a full upstream clone if you need to rebuild `fsm.exe` too,
  which the vendor snapshot does not support — see the vendor README).
- `patches/` — for-reference copies of the upstream C++ patches
  required by FSM_Julia. They are **pre-applied** in
  `../vendor/Barnes2020-FillSpillMerge/`, so you only need to interact
  with them if you're refreshing the vendor from a clean upstream clone
  or auditing the diff vs. upstream. See `patches/README.md`.

## Regenerating oracle data

```sh
# 1. Build the helper binaries against the vendored C++ snapshot.
./build.sh

# 2. Regenerate the Phase 2 dephier oracle for each test case:
for d in ../test/test_cases/case_*; do
  ./dh_dump.exe "$d/input.tif" 0.0 "$d/expected-dh.txt"
done

# 3. Regenerate the Phase 4 random terrain batch (overwrites the
#    existing 70 .tif files under test/test_cases/random/):
./gen_random_terrains.exe ../test/test_cases/random
```

Regenerating `expected-wtd.tif` is not possible against the vendored
snapshot alone (that file is produced by upstream's `fsm.exe`, which
needs the full upstream source tree). To refresh those, clone the full
upstream repo, apply the two patches in `patches/`, build via CMake,
and run `fsm.exe` per case.
