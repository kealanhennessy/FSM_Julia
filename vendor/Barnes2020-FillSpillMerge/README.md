# Vendored Barnes2020-FillSpillMerge (minimal subset)

This directory is a **minimal vendored snapshot** of the upstream
[Barnes2020-FillSpillMerge](https://github.com/r-barnes/Barnes2020-FillSpillMerge)
C++ project (plus its two transitive submodules
[dephier](https://github.com/r-barnes/Barnes2019-DepressionHierarchy)
and [richdem](https://github.com/r-barnes/richdem)), included so that
this branch of FSM_Julia is self-contained: cloning this repo gives you
everything you need to rebuild the test-data tools and re-run the
oracle comparisons against the reference C++ algorithm, without
fetching anything else.

## Scope

Only what `FSM_Julia/tools/` needs to build is here. Specifically:

| Path | Why |
|---|---|
| `include/fsm/fill_spill_merge.hpp` | The C++ FSM algorithm header (the file `src/fill_spill_merge.jl` was ported from). Not directly consumed by the test tools but kept as the canonical reference. |
| `submodules/dephier/include/dephier/*.hpp` | dephier headers consumed by `tools/dh_dump.cpp`. |
| `submodules/dephier/submodules/richdem/include/**` | richdem headers transitively included by dephier and the test tools. Whole `include/` tree is here (~2.4 MB) so the include graph never breaks on a refactor of the upstream's internal headers. |
| `submodules/dephier/submodules/richdem/src/terrain_generation/{PerlinNoise.*,terrain_generation.cpp}` | richdem's perlin implementation, linked into `tools/gen_random_terrains.exe`. |

What is **not** here:

- The upstream `src/main.cpp`, so you cannot rebuild `fsm.exe`
  end-to-end. The committed `test/test_cases/case_*/expected-wtd.tif`
  files were produced by `fsm.exe` against the patched upstream and
  remain the oracle for the Julia port's bit-exact comparison.
- The upstream `unittests/`, `paper/`, scaling tests, doc trees, CMake
  files, and richdem wrappers / Python bindings — none of which are
  required for the Julia tests to run.
- Submodule `.git` directories. The vendored snapshot is flat; you
  cannot run `git log` inside `submodules/dephier/`.

If you want to regenerate `expected-wtd.tif` (the only thing the
vendored copy doesn't support), pull the full upstream + apply the
patches in `FSM_Julia/tools/patches/` and run its CMake build.

## Source

The vendored files are taken from these exact upstream commits:

| Repo | Commit |
|---|---|
| Barnes2020-FillSpillMerge | `1c499ea475c09b9f4c5da74ee5cc995de169db63` |
| Barnes2019-DepressionHierarchy (dephier) | `411f7d4ad344d74447b47cf9eb85acd536a4d8a1` |
| richdem | `415032db2f30372111e4cfd37f046e7542ed66f3` |

## Patches applied

Two patches are **pre-applied** in this vendored copy. They are reproduced
verbatim in `FSM_Julia/tools/patches/` as standalone files for audit
(diff vs. upstream) and so anyone refreshing this vendor from a clean
upstream clone can re-derive the same state.

1. **`dephier-deterministic-outlets.patch`** — extends the outlet-sort
   comparator in
   `submodules/dephier/include/dephier/dephier.hpp` to break ties on
   `(depa, depb, out_cell)` after the elevation key. Without this,
   `std::sort`'s instability + `unordered_map` iteration order make the
   resulting meta-depression tree implementation-defined, which blocks
   bit-exact cross-language comparison. Physically the algorithm is
   insensitive to which equally-low outlet is picked.

2. **`richdem-gdal-cslconstlist.patch`** — changes the `ProcessMetadata`
   parameter type in
   `submodules/dephier/submodules/richdem/include/richdem/common/Array2D.hpp`
   from `char**` to `CSLConstList`, matching the signature
   GDAL ≥ 3.7 exposes for `GDALDataset::GetMetadata()`. The upstream
   compiles fine against older GDAL but not against the version on the
   port author's machine.

Neither patch was pushed upstream — see `FSM_Julia/tools/patches/`
header notes and the project memory for context.

## Licenses

Each upstream repo's LICENSE file is preserved at the corresponding path:

- `LICENSE` — MIT, Copyright (c) 2020 Richard Barnes (FSM)
- `submodules/dephier/LICENSE` — MIT, Copyright (c) 2020 Richard Barnes (dephier)
- `submodules/dephier/submodules/richdem/LICENSE.txt` — **GPL v3** (richdem)

richdem being GPL3 means the vendored subdirectory at
`submodules/dephier/submodules/richdem/` carries GPL3 obligations
(preserve LICENSE.txt, offer source on redistribution). The Julia code
under `FSM_Julia/src/` doesn't link against richdem at runtime — it
only reads data files produced by separately-invoked C++ binaries
(`tools/dh_dump.exe`, `tools/gen_random_terrains.exe`) — so the
"separate program, mere aggregation" reading of the GPL applies. The
Julia port's own licensing is up to FSM_Julia's top-level LICENSE (if
one is set).
