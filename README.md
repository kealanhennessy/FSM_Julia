# FSM_Julia

A Julia port of the C++ Fill-Spill-Merge (FSM) algorithm of
Barnes ([2020](https://github.com/r-barnes/Barnes2020-FillSpillMerge)).
The port is line-by-line faithful to the upstream so that outputs match
the reference `fsm.exe` bit-for-bit on the committed test cases.

Running the test suite is **fully self-contained** — clone the repo,
have Julia, done. No C++ toolchain, no external repositories. All the
reference data the tests compare against is committed here.

Regenerating that reference data (a rare maintenance task) is the only
thing that needs the upstream C++. This repo does **not** bundle a copy
of it; instead it ships the exact patches and step-by-step instructions
for setting up your own upstream clone. See Testing → Option B.

## Quick start

```sh
git clone <this repo>
cd FSM_Julia
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Expected: `Pkg.test()` reports `1700 / 1700 pass` in roughly 7 seconds.
Nothing besides Julia is needed for this — see Testing → Option A.

---

# Testing

There are two ways to exercise this repository, and they are completely
independent of each other:

- **Option A — Run the test suite (`Pkg.test()`).** The normal path.
  Pure Julia. Compares the port against reference data files that are
  already committed to the repo. **No C++ compiler, no GDAL, nothing
  but Julia.** This is what you almost always want.

- **Option B — Regenerate the oracle data.** Only needed if you have
  changed the algorithm or a test fixture and need to rebuild the
  reference data the suite compares against. This is the only path that
  compiles C++.

The relationship between them: Option B *produces* the frozen reference
files (`expected-dh.txt`, the random terrain `.tif`s, the Priority-Flood
oracles, and `expected-wtd.tif`); Option A *consumes* them. The `tools/`
programs and the patches in `tools/patches/` exist solely to support
Option B, against an upstream clone you set up yourself. Nothing in
Option B is touched by Option A.

## Option A — Run the test suite

### A.1 Prerequisites

- **Julia ≥ 1.11** (the reference machine runs 1.12.6).
- Network access on first run so `Pkg.instantiate()` can fetch the
  Julia dependencies (`ArchGDAL`, `ArgParse`, `DataStructures`, `Test`)
  — about 140 transitive packages, a few hundred MB, cached in the
  Julia depot afterwards.

You do **not** need a C++ toolchain, GDAL, CMake, or anything else. All
reference data is committed under `test/test_cases/` (~3.5 MB: five
hand-crafted oracle cases plus seventy generated terrains).

### A.2 Steps

```sh
git clone <repo-url> FSM_Julia
cd FSM_Julia

# Resolve Project.toml/Manifest.toml to the exact pinned dependency
# versions (one-time; a minute or two on first run):
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run everything:
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected final output:

```
Test Summary:  | Pass  Total  Time
FillSpillMerge | 1700   1700  ~7s
```

### A.3 What it actually does

`Pkg.test()` runs `test/runtests.jl`, which loads the Julia port and,
for each test group below, either checks an invariant directly or reads
a committed reference file and diffs the port's output against it. The
reference files are static inputs to the test run; nothing regenerates
them and no C++ is involved.

### A.4 Test groups

| Group | Tests | What it checks |
|---|---|---|
| **Core utilities** | ~30 | D8 offset/constant consistency, the `fp_eq`/`fp_le`/`fp_ge` floating-point comparison helpers, and `BucketFillFromEdges` (ocean flood-fill). |
| **Data structures** | ~30 | `DisjointDenseIntSet` (union/find with the merge-into-parent variant), the LIFO-min priority queue, and `Depression`/`Outlet` struct defaults. |
| **Depression-hierarchy oracles** | ~250 | `get_depression_hierarchy` run on each of the 5 hand-crafted cases, compared **bit-exact** against `expected-dh.txt` (label grid, flow directions, and every field of every `Depression`). |
| **Water-table-depth oracles** | 15 | Full FSM run on each of the 5 cases, compared against `expected-wtd.tif` produced by the reference C++ `fsm.exe`. Tolerance 1e-9 for cases 1–4 (exact rationals), 1e-6 for the 100×100 case. |
| **Ported C++ unit tests** | ~32 | Direct 1:1 ports of every `TEST_CASE` in the upstream `unittests/fsm_tests.cpp` (DepressionVolume, DetermineWaterLevel, MoveWaterIntoPits, BackfillDepression, FillDepressions and its sub-cases, the two "PQ Issue" regressions). Deterministic; no random input. |
| **Randomized property tests** | ~1130 | Five properties iterated over the 70 generated terrains: surface water is conserved; FSM is idempotent under re-application; incremental vs. bulk water input agree; `MoveWaterIntoPits` is consistent on re-seeding; and a heavy-flood run matches a frozen, trusted-C++ Priority-Flood (Zhou 2016) oracle. |
| **Edge cases** | ~70 | 1×1 land (throws — needs ocean), all-ocean (throws), lone interior cell absorbed into ocean, monotonic slope (no depressions), a pit filled exactly to its brim, and a full run on `Float32` elevations. |

(The group sizes are approximate because several groups loop over the
70 terrains or per-depression fields, multiplying the assertion count.)

## Option B — Regenerate the oracle data

You only need this if you changed the algorithm or a test fixture and
must rebuild the reference data Option A compares against. It is the
only path that compiles C++, and the only thing that needs the upstream
C++ project. This repo intentionally does **not** bundle the upstream;
you set up your own clone once (B.2), and from then on regenerating any
oracle is mechanical.

### B.1 What the reference data is, and what regenerates it

Every file below is committed; this is just what *produces* each when
you regenerate. All four require the upstream-clone setup in B.2.

| Reference file(s) | Produced by |
|---|---|
| `test/test_cases/case_*/expected-dh.txt` | `tools/dh_dump.exe` |
| `test/test_cases/random/*.tif` (70 files) | `tools/gen_random_terrains.exe` |
| `test/test_cases/random_pf/*.tif` (35 files) | `tools/pf_dump.exe` |
| `test/test_cases/case_*/expected-wtd.tif` | upstream's own `fsm.exe` |

The first three come from the small helper programs in `tools/`; the
last from the upstream's own `fsm.exe`. All four are built from the
same patched upstream clone.

### B.2 One-time setup: clone, pin, patch the upstream

You need a C++17 compiler, GDAL (with `gdal-config` on `PATH`), CMake,
and `git`. On macOS: `brew install gdal cmake`. On Debian/Ubuntu:
`apt install g++ libgdal-dev cmake git`.

The patches were generated against specific upstream commits. Check the
parent repo out at that commit so its submodules resolve to the matching
pinned revisions, then apply the two patches:

```sh
# Pinned upstream commits (parent repo pins the two submodules):
#   Barnes2020-FillSpillMerge          1c499ea475c09b9f4c5da74ee5cc995de169db63
#   Barnes2019-DepressionHierarchy     411f7d4ad344d74447b47cf9eb85acd536a4d8a1  (submodules/dephier)
#   richdem                            415032db2f30372111e4cfd37f046e7542ed66f3  (…/richdem)

git clone https://github.com/r-barnes/Barnes2020-FillSpillMerge.git
cd Barnes2020-FillSpillMerge
git checkout 1c499ea475c09b9f4c5da74ee5cc995de169db63
git submodule update --init --recursive   # pulls the pinned submodule commits

# Apply this port's two required patches (rationale in "C++ patches",
# below). $FSM_JULIA = path to your clone of THIS repo.
git -C submodules/dephier apply \
  "$FSM_JULIA/tools/patches/dephier-deterministic-outlets.patch"
git -C submodules/dephier/submodules/richdem apply \
  "$FSM_JULIA/tools/patches/richdem-gdal-cslconstlist.patch"

export FSM_CPP="$PWD"          # used by tools/build.sh in B.3
```

Newer upstream HEAD may have moved the patched lines; checking out the
pinned commit above is what guarantees the patches apply cleanly. How
you make GDAL/CMake happy on your platform (package versions,
cross-arch, etc.) is left to you — that's deliberately *your* setup, not
something this repo tries to abstract.

### B.3 Build everything against that clone

```sh
# (a) The three helper programs. FSM_CPP must point at the patched
#     clone from B.2. tools/build.sh uses `c++` from PATH and the
#     host's native architecture; no machine-specific paths are baked
#     in. Two optional overrides for unusual setups:
#       CXX=    a specific compiler (e.g. /usr/local/opt/llvm/bin/clang++)
#       ARCH=   cross-target, needed ONLY if your GDAL's arch differs
#               from the host (classic case: Apple Silicon with an
#               x86_64 Homebrew GDAL -> ARCH=x86_64)
cd "$FSM_JULIA"
tools/build.sh                 # builds tools/{dh_dump,pf_dump,gen_random_terrains}.exe

# (b) The upstream's own fsm.exe, via its CMake (for expected-wtd.tif).
#     Match the same arch as your GDAL if cross-building.
cmake -S "$FSM_CPP" -B "$FSM_CPP/build" -DUSE_GDAL=ON
cmake --build "$FSM_CPP/build" --target fsm.exe
```

What the helper programs do:

- **`tools/dh_dump.exe`** — runs the C++ ocean-labeling +
  `GetDepressionHierarchy` on a GeoTIFF and writes the label grid,
  flow-directions grid, and full `Depression` vector to the plain-text
  format the depression-hierarchy oracle tests parse.
- **`tools/pf_dump.exe`** — runs the C++ Zhou (2016) Priority-Flood on
  a terrain GeoTIFF and writes the filled DEM back out as a Float64
  GeoTIFF. The trusted reference the heavy-flooding property test
  compares FSM against.
- **`tools/gen_random_terrains.exe`** — uses richdem's `perlin` to emit
  the deterministic batch of 70 small/large × integer/float terrains.
  Bumping `MASTER_SEED` in the source yields a fresh batch.

### B.4 Regenerate the oracles

```sh
cd "$FSM_JULIA"

# Depression-hierarchy oracles (one expected-dh.txt per case):
for d in test/test_cases/case_*; do
  tools/dh_dump.exe "$d/input.tif" 0.0 "$d/expected-dh.txt"
done

# The 70-terrain random batch (overwrites test/test_cases/random/*.tif):
tools/gen_random_terrains.exe test/test_cases/random

# The Priority-Flood oracles, one per float terrain (35 files). Run
# this whenever the random terrains are regenerated:
mkdir -p test/test_cases/random_pf
for f in test/test_cases/random/small_float_*.tif \
         test/test_cases/random/large_float_*.tif; do
  tools/pf_dump.exe "$f" "test/test_cases/random_pf/$(basename "$f")"
done

# The water-table-depth oracles. run_all.sh (and each case's run.sh)
# requires FSM_EXE — there is no default path. Each case's run.sh
# already encodes that case's ocean_level / --swl (see Test cases).
FSM_EXE="$FSM_CPP/build/fsm.exe" test/test_cases/run_all.sh
```

### B.5 Re-verify

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Note: regenerating the random terrains rewrites the `.tif` files; their
*elevation data* is deterministic from `MASTER_SEED`, but GeoTIFF
headers embed a timestamp, so the files will show as changed in git
even when the numbers are identical. The same applies to the
Priority-Flood oracle `.tif`s.

---

# Test cases

The five hand-crafted cases under `test/test_cases/case_*/` each
contain:

- `input.asc` — Arc/Info ASCII grid of the DEM (human-readable source).
- `input.tif` — same DEM as a Float64 GeoTIFF (the format FSM consumes).
- `expected-wtd.tif` — water-table depth from the reference C++ binary.
- `expected-hydrologic-surface-height.tif` — `wtd + topo`.
- `expected-dh.txt` — the dumped depression hierarchy (regenerable; see
  Option B).
- `run.sh` — reproduces the case from `input.asc` via the upstream
  binary. (Case 5 also has `generate.py` to re-derive its input.)

### Algorithm contract exercised by the cases

The reference binary is invoked as
`fsm.exe <topo.tif> <out_prefix> <ocean_level> --swl <X>`:

- **`ocean_level`** — FSM bucket-fills inward from the grid edges,
  marking edge-connected cells with elevation ≤ `ocean_level` as OCEAN.
  Any NoData cell also becomes OCEAN regardless of position. FSM errors
  if there are no OCEAN cells, so every case has at least one NoData
  cell.
- **`--swl X`** — constant surface-water level: every non-ocean cell
  starts with X units of standing water.
- **Output** — `<prefix>-wtd.tif` (water-table depth) and
  `<prefix>-hydrologic-surface-height.tif` (`wtd + topo`). Ocean cells
  have `wtd = 0`.

### The cases

| # | Size | `ocean_level` | `--swl` | Tol | Purpose |
|---|---|---|---|---|---|
| 1 | 5×5 | −100 | 1.0 | 1e-9 | Trivial trough; one NoData ocean cell. Sanity check + corner-ocean drainage. |
| 2 | 10×10 | 0 | 1.0 | 1e-9 | Single bowl-shaped depression, NoData border ocean. No overflow. |
| 3 | 10×10 | 0 | 1.5 | 1e-9 | Two adjacent depressions sharing a sill. Forces overflow + merge. |
| 4 | 10×10 | 0 | 1.0 | 1e-9 | Three-level hierarchy: two pits in a bowl in a plateau. |
| 5 | 100×100 | 0 | 0.5 | 1e-6 | Synthetic 1/f-noise terrain; many depressions; full-scale exercise. |

The 70 files under `test/test_cases/random/` are perlin terrains
(30 small + 5 large, each in integer-truncated and float variants) used
only by the randomized property tests. Four of the five property tests
assert invariants and need no committed "expected" output. The fifth —
heavy-flooding-vs-Priority-Flood — does compare against a frozen
reference: `test/test_cases/random_pf/` holds one trusted-C++
Priority-Flood–filled DEM per *float* terrain (35 files, same
basenames), produced by `tools/pf_dump.exe`.

---

# Testing-support code

None of the code below is part of the Fill-Spill-Merge algorithm. It
was written purely to make the port testable: parsing reference dumps,
bridging C++ and Julia indexing conventions, and loading fixtures. It
lives in `test/` (Julia harness helpers) and `tools/` (C++ programs
that produce the frozen reference data). Each item below gives the
signature, what it does, and *why* it had to exist.

## Harness helpers in `test/runtests.jl`

These are defined in the test entry point itself because they are the
glue between the suite and the port.

- **`read_tif(path) -> Matrix`** — Reads band 1 of a GeoTIFF via
  ArchGDAL and, if the band declares a non-NaN NoData value, replaces
  those cells with `NaN`. *Why:* input topographies encode the ocean
  border as a NoData value; the algorithm and the rest of the harness
  expect ocean as `NaN`. For the `fsm.exe` reference outputs this is a
  deliberate no-op — the C++ binary writes a numeric `0` at ocean
  cells (NoData is some other sentinel), so the oracle is read back
  faithfully rather than having its zeros turned into `NaN`.

- **`compute_julia_wtd(topo, ocean_level, swl) -> Matrix`** — The
  single integration seam between the harness and the port. It mirrors
  the setup in the upstream `src/main.cpp`: build the label/flowdirs
  grids and ocean labelling via `prepare_label_and_flowdirs`, run
  `get_depression_hierarchy`, initialise `wtd` to `swl` on land and `0`
  on ocean cells, then run `fill_spill_merge!`. Returns the resulting
  water-table-depth grid. *Why:* the water-table-depth oracle tests
  need exactly the same pre-processing the C++ `main.cpp` does around
  the algorithm, isolated in one place so the test body is just
  "compute vs. expected".

- **`CASES`, `TEST_CASES_DIR`, `PORT_READY`** — `CASES` is the table of
  the five hand-crafted cases and their `ocean_level` / `swl` /
  tolerance (the parameters each case's reference output was generated
  with). `TEST_CASES_DIR` is the fixtures path. `PORT_READY` is a flag
  that gates the water-table-depth oracle tests: when `true` (its
  current, permanent value) those tests run for real; when `false` they
  become `@test_skip`. It was the switch used while the port was being
  built and is left in place as an explicit, greppable marker of that
  gate.

## Oracle parsing + comparison in `test/dephier_oracle.jl`

This file exists entirely to consume the depression-hierarchy dump
produced by `tools/dh_dump.exe` and compare it to the Julia port.

- **`read_dh_dump(path) -> OracleDH`** — Parses the plain-text dump
  format (`WIDTH` / `HEIGHT` / `NDEP` header, then the `LABEL` grid,
  the `FLOWDIRS` grid, then one line per `Depression`). Returns an
  `OracleDH`. *Why:* the C++ side can't hand Julia a struct, so the
  hierarchy is serialised to text by `dh_dump.cpp` and reconstructed
  here.

- **`OracleDH`** — Immutable struct
  `{W::Int, H::Int, label::Matrix{dh_label_t}, flowdirs::Matrix{Int8},
  deps::Vector{Depression{Float64}}}`. The parsed reference hierarchy in
  the same field types the port produces, so comparison is direct.

- **`prepare_label_and_flowdirs(topo_in, ocean_level)
  -> (label, flowdirs, topo_clean)`** — Reproduces the ocean-labelling
  preamble of the upstream `main.cpp`: copies `topo_in`, replaces every
  `NaN` with `-Inf`, allocates `label` (all `NO_DEP`) and `flowdirs`
  (all `NO_FLOW`), bucket-fills ocean inward from the edges at
  `ocean_level`, then marks every `-Inf` cell (and anything the bucket
  fill already flagged) as `OCEAN`. *Why the `NaN → -Inf` substitution:*
  the C++ `Array2D` keeps the raw NoData numeric value (e.g. `-9999`),
  which orders correctly in the dephier priority queue. `read_tif`
  represents ocean as `NaN`, and `NaN` breaks heap ordering (every
  comparison is false), so ocean cells would not pop first. `-Inf` is
  order-equivalent to a very negative NoData and keeps the OCEAN
  wavefront popping ahead of land — without it the hierarchy would
  differ from the C++. This helper is used both by `compute_julia_wtd`
  and by every edge-case test.

- **`diff_grid(name, jl, c) -> (matched::Bool, count::Int, first)`** —
  Element-by-element comparison of two grids using `isequal`, so that
  `NaN` compares equal to `NaN` (ordinary `==` would make every ocean
  cell a spurious mismatch). Returns whether they matched, the
  mismatch count, and a named tuple describing the first mismatch (or
  `nothing`). *Why:* used by the depression-hierarchy oracle tests to
  diff the `label` and `flowdirs` grids and report the first divergence
  precisely.

- **`diff_depression(jl, c) -> Vector{Tuple{Symbol,Any,Any}}`** —
  Compares all sixteen `Depression` fields with `isequal` and returns
  `(field, julia_value, cpp_value)` for each mismatch. *Why:* lets the
  oracle test pinpoint exactly which field of which depression diverged
  rather than just "structs differ".

- **`_parse_float(s) -> Float64`** — Parses a float token, special-casing
  the `nan` / `-nan` / `inf` / `+inf` / `-inf` spellings a C++ ostream
  emits, otherwise `parse(Float64, s)`. *Why:* Julia's
  `parse(Float64, ...)` handling of those spellings has varied across
  versions; the dump contains `inf`/`nan` for sentinel elevations, so
  parsing has to be explicit to stay version-independent.

- **`_to_julia_flat(i) -> UInt32`** — `i == NO_VALUE ? NO_VALUE : i + 1`.
  Converts a 0-based C++ flat cell index to a 1-based Julia linear
  index, leaving the `NO_VALUE` sentinel untouched. *Why:* this is the
  *only* place the 0-based↔1-based bridge happens. The algorithm port
  is 1-based throughout; the conversion is confined to the oracle
  boundary so it can't leak into the algorithm.

## Fixture helpers in `test/fsm_unit_tests.jl`

Used by the ported deterministic unit tests and reused by the edge
tests.

- **`visual_grid(m) = permutedims(m)`** — Transposes a matrix literal.
  *Why:* the upstream C++ tests write grid fixtures visually (each
  source line is a row of the grid). Julia's `M[i, j]` is row `i`, col
  `j`, but the project convention is `M[x, y]` = column `x`, row `y`
  (to mirror C++ `Array2D`'s `M(x, y)`). Writing the fixture in visual
  form and wrapping it in `visual_grid` yields a `(W, H)` array whose
  `topo[x, y]` equals the C++ source's `topo(x, y)` — so fixtures can
  be transcribed from the upstream verbatim without hand-transposing.

- **`xy_to_flat(W, H, x, y) -> Int`** — `LinearIndices((W, H))[x+1, y+1]`.
  Converts a C++ 0-based `(x, y)` coordinate (as written in the
  upstream test source, e.g. `topo.xyToI(4, 2)`) to a 1-based Julia
  flat linear index. *Why:* pit/outlet cells in the ported fixtures
  are specified with the upstream's 0-based coordinates; this keeps the
  transcription mechanical and the off-by-one in one audited place.

- **`set_edges_ocean!(label)`** — Sets the outer one-cell ring of
  `label` to `OCEAN`, in place. Mirrors the upstream's
  `Array2D::setEdges(OCEAN)`. *Why:* several upstream tests (and the
  random-terrain setup) establish the ocean by labelling the border
  ring directly rather than by elevation; this reproduces that exactly.

## Fixture helpers in `test/fsm_random_tests.jl`

Drive the property tests over the 70 generated terrains.

- **`RANDOM_TERRAINS_DIR`** — Path constant for
  `test/test_cases/random/`.

- **`_list_terrains(pattern) -> Vector{String}`** — Returns the sorted
  list of terrain paths whose filename starts with `pattern` (e.g.
  `"small_int_"`, `"large_float_"`); empty vector if the directory is
  absent. *Why:* lets each property test select the relevant subset
  (integer vs. float, small vs. large) the corresponding upstream
  `TEST_CASE` used.

- **`_load_random_terrain(path) -> Matrix{Float64}`** — Loads band 1 as
  a `Float64` matrix **without** the NoData→`NaN` substitution
  `read_tif` does. *Why:* the generated terrains encode the ocean ring
  as a numeric `-1.0` (the generator's declared NoData is an unused
  `-9999`). The property tests need that `-1.0` ring preserved
  numerically — converting it to `NaN` here would change the input the
  algorithm sees relative to how the batch was generated.

- **`_init_label_flowdirs(W, H) -> (label, flowdirs)`** — Allocates
  `label` (all `NO_DEP`) and `flowdirs` (all `NO_FLOW`) and applies
  `set_edges_ocean!`. *Why:* every property test starts from the same
  border-ocean setup the upstream randomized tests use; this removes
  that boilerplate from each test body.

## C++ reference-data programs in `tools/`

Fully described under Testing → Option B (B.3). Summarised here as
testing-support code. All three are compiled only by `tools/build.sh`
(Option B) and never touched by `Pkg.test()` (Option A) — they only
*produce* the frozen reference data Option A consumes.

- **`tools/dh_dump.cpp`** — links your upstream C++ dephier clone (B.2),
  runs the same ocean-labelling + `GetDepressionHierarchy` the upstream
  does, and serialises the label grid, flow-directions grid, and every
  `Depression` to the text format `read_dh_dump` parses. Producer of
  the depression-hierarchy oracle.

- **`tools/pf_dump.cpp`** — links your upstream C++ richdem clone's
  Zhou (2016) Priority-Flood, runs it on a terrain GeoTIFF, and writes
  the filled DEM back out as a Float64 GeoTIFF. Producer of the
  Priority-Flood
  oracle under `test/test_cases/random_pf/`. *Why this exists at all:*
  the heavy-flooding property test needs an *independent* algorithm to
  cross-check FSM against (two unrelated algorithms agreeing is far
  stronger evidence than FSM agreeing with itself). Earlier this
  reference was a Julia re-implementation of Priority-Flood, but a
  same-author re-implementation has its own unaudited failure modes and
  no oracle of its own. Dumping the *trusted upstream C++* Priority-Flood
  output as a frozen oracle — exactly how the dephier and
  water-table-depth references work — keeps the cross-algorithm
  independence while grounding the reference in code that is already
  trusted. (The Julia Priority-Flood port was consequently removed.)

- **`tools/gen_random_terrains.cpp`** — links your upstream C++ richdem
  clone's `perlin` and emits the deterministic 70-terrain batch the
  property tests iterate over.

---

# Reference system

The only configuration the port has been exercised on end-to-end.
Anything close should work; the exact numbers are recorded so
reproducibility issues can be triaged against a known-good baseline.

| Component | Version |
|---|---|
| Hardware | Apple M1 Pro |
| RAM | 16 GB |
| OS | macOS Tahoe 26.3.1 |
| Julia | 1.12.6 |
| C++ compiler (Option B only) | Apple clang (`c++`), or Homebrew LLVM clang 22.1.5 |
| GDAL (Option B only) | 3.13.0 (Homebrew, x86_64) |

Only Option B touches a C++ toolchain at all; Option A is pure Julia.
Nothing in the repo hardcodes a path on this machine. `tools/build.sh`
uses `c++` from `PATH` and the host's native architecture by default.

This particular machine has one wrinkle worth recording: its Homebrew
(and therefore its GDAL) lives under `/usr/local` — the Intel prefix,
running via Rosetta on Apple Silicon — so GDAL is an x86_64 build while
the host is arm64. On a setup like that, the helper tools must be
cross-built for x86_64:

```sh
ARCH=x86_64 tools/build.sh
```

On a machine where GDAL matches the host architecture (a plain
`brew install gdal` on Apple Silicon, or any Linux box), no override is
needed — just `tools/build.sh`. The Julia tests (Option A) don't care
about the helper binaries' architecture; they only consume the data
those binaries produced.

---

# Repository layout

```
src/                          The Julia port (the algorithm itself)
  FillSpillMerge.jl           Top-level module; includes the rest
  constants.jl                D8 offsets, OCEAN / NO_DEP sentinels
  fp_compare.jl               fp_eq / fp_le / fp_ge (1e-6 tolerance)
  bucket_fill.jl              BucketFillFromEdges
  types.jl                    Depression / Outlet structs
  disjoint_set.jl             DisjointDenseIntSet
  priority_queue.jl           LIFOMinPriorityQueue
  dephier.jl                  GetDepressionHierarchy
  fill_spill_merge.jl         FillSpillMerge driver + 9 helpers

test/
  runtests.jl                 Test runner entry point
  dephier_oracle.jl           expected-dh.txt parser + shared helpers
  fsm_unit_tests.jl           Ported deterministic C++ unit tests
  fsm_random_tests.jl         Randomized property tests
  fsm_edge_tests.jl           Edge cases + Float32
  test_cases/
    case_01_trough/ … case_05_perlin_100x100/   Hand-crafted oracle cases
    random/                                       70 perlin terrains
    random_pf/                                    35 C++ Priority-Flood oracles

tools/                        C++ helpers — only used by Option B
  dh_dump.cpp                 Dumps GetDepressionHierarchy to text
  pf_dump.cpp                 Dumps a C++ Zhou2016 Priority-Flood fill
  gen_random_terrains.cpp     Emits the 70-terrain batch
  build.sh                    Builds all three against $FSM_CPP (your
                              upstream clone — see Testing → Option B)
  patches/                    The two required upstream patches + the
                              pinned commits they apply against
```

No upstream C++ is bundled in this repo. The `tools/` programs and
upstream `fsm.exe` build against a clone you set up yourself (Testing →
Option B); the patches and exact commit SHAs needed to reproduce that
setup live in `tools/patches/`.

---

# Upstream C++ provenance

The reference data (and anything you regenerate in Option B) is produced
from these exact upstream commits — the parent repo pins its two
submodules, so checking it out at the commit below resolves the matching
dephier/richdem revisions:

| Repo | Commit |
|---|---|
| [Barnes2020-FillSpillMerge](https://github.com/r-barnes/Barnes2020-FillSpillMerge) | `1c499ea475c09b9f4c5da74ee5cc995de169db63` |
| [Barnes2019-DepressionHierarchy](https://github.com/r-barnes/Barnes2019-DepressionHierarchy) (`submodules/dephier`) | `411f7d4ad344d74447b47cf9eb85acd536a4d8a1` |
| [richdem](https://github.com/r-barnes/richdem) (`…/richdem`) | `415032db2f30372111e4cfd37f046e7542ed66f3` |

Two small patches (in `tools/patches/`) must be applied to that clone
before the helper tools or `fsm.exe` will build/behave correctly —
walked through in Testing → Option B (B.2) and justified below. Nothing
from these repos is redistributed here; only the unified-diff patches we
authored are committed.

---

# C++ patches: what they change and why they were necessary

Two small modifications to the upstream C++ are required before it can
serve as a bit-exact oracle for the Julia port. They are committed as
standalone unified diffs in `tools/patches/` and applied to your own
upstream clone during Option B setup (B.2); neither was pushed upstream
(the upstream repo belongs to a different author). The justification for
each follows so the modifications can be defended on review.

### Patch 1 — deterministic outlet tie-break in dephier

**File:** `submodules/dephier/include/dephier/dephier.hpp`, in
`GetDepressionHierarchy`, the outlet sort (~line 620).

**Upstream code** sorts candidate outlets by elevation only:

```cpp
std::sort(outlets.begin(), outlets.end(),
  [](const Outlet<elev_t> &a, const Outlet<elev_t> &b){
    return a.out_elev < b.out_elev;
  });
```

**Patched code** adds a total-order tie-break:

```cpp
std::sort(outlets.begin(), outlets.end(),
  [](const Outlet<elev_t> &a, const Outlet<elev_t> &b){
    return std::tie(a.out_elev, a.depa, a.depb, a.out_cell)
         < std::tie(b.out_elev, b.depa, b.depb, b.out_cell);
  });
```

**Why it is necessary.** When two outlets sit at *exactly* the same
elevation, the upstream comparator treats them as equivalent, so their
final relative order is decided by two implementation-defined factors:

1. `std::sort` is **not stable** — equal elements may be reordered
   arbitrarily by the standard library's sort implementation.
2. The outlets are copied out of a `std::unordered_map`, whose
   **iteration order is unspecified** by the C++ standard.

A single compiled binary is self-consistent (same result every run on
that build), so the upstream never needed to pin this down. But a
re-implementation in another language cannot reproduce libc++/libstdc++
sort internals or a particular `unordered_map` bucket layout, so the
two implementations diverge on equal-elevation ties.

**Why the divergence matters for testing but not for the science.**
Different tie orderings produce *different binary trees that encode the
same physical depression structure*:

- Leaf depressions (one per pit cell) are unaffected — the watershed
  traversal that creates them is fully deterministic.
- Meta-depressions are formed by merging adjacent depressions when
  their shared outlet is processed. With several equally-low outlets, a
  different processing order merges a different pair first, yielding a
  tree of identical physical content but different shape.
- `CalculateMarginalVolumes` walks each cell up its parent chain until
  it finds a depression whose `out_elev` exceeds the cell's elevation.
  A different tree shape changes which meta-depression a cell is
  attributed to, so per-depression `cell_count` / `total_elevation` /
  `dep_vol` can differ between orderings — even though the totals over
  the whole tree are identical.

So the patch does not "fix a bug": Barnes (2020) does not specify a
tie-break because, scientifically, there is nothing to specify — all
orderings yield equivalent hierarchies with identical whole-tree
results. The patch only removes an implementation-defined degree of
freedom so that two independent implementations can be compared
**bit-for-bit** (the depression-hierarchy oracle). The Julia port uses
the matching key `sort!(outlets, by = o -> (o.out_elev, o.depa,
o.depb, o.out_cell))`; with both sides sharing the comparator the
hierarchies are identical.

**Cost / risk.** One `std::tie` of three already-present integer fields
per sort comparison, dominated by the sort itself — negligible runtime
cost, zero change to any scientifically-meaningful output.

### Patch 2 — GDAL `CSLConstList` compatibility in richdem

**File:**
`submodules/dephier/submodules/richdem/include/richdem/common/Array2D.hpp`,
the `ProcessMetadata` free function (~line 42).

**Change:** the parameter type goes from `char **metadata` to
`CSLConstList metadata` (one token).

**Why it is necessary.** richdem calls `ProcessMetadata` with the
result of `GDALDataset::GetMetadata()`. In a current GDAL that method
returns `CSLConstList` (an alias for `const char* const*`), not the
`char**` the upstream signature expects. This is directly verifiable in
the GDAL installed on the reference machine — `gdal_majorobject.h:79`
in GDAL 3.13.0 declares:

```cpp
virtual CSLConstList GetMetadata(const char *pszDomain = "");
```

Passing a `const char* const*` to a `char**` parameter is a
const-correctness violation and a hard compile error: the upstream
simply **does not compile** against this GDAL, so there is no way to
build `dh_dump.exe` / `gen_random_terrains.exe` (or upstream's own
`fsm.exe`) without the change. (The `char**` signature was valid
against the older GDAL the upstream was originally written for; the
return type was narrowed to `CSLConstList` somewhere in the GDAL 3.x
series — the exact minor version isn't pinned here because the
justification rests only on the GDAL actually in use, 3.13.0.)

**Why it is safe.** This is a pure build-compatibility fix. The
function body is unchanged; `CSLConstList` is exactly the type GDAL now
hands in, so the patched signature is *more* correct than the original,
not a workaround. There is no behavioural change — it does not touch
any algorithm path, only lets the existing code compile against a
modern GDAL. (It was also independently needed early on just to get the
reference `fsm.exe` to build at all.)

---

# Licenses

**No third-party source is redistributed in this repository.** The
upstream C++ is not bundled — it is cloned by the user during Option B.
That deliberately keeps the licensing simple:

- The Julia port (`src/`, `test/`) and the helper-tool sources
  (`tools/*.cpp`, `tools/build.sh`) are this project's own work; their
  licensing is left to a top-level `LICENSE` file (not yet set on this
  branch).
- `tools/patches/*.patch` are short unified diffs *we* authored. One
  targets dephier (MIT, Copyright (c) 2020 Richard Barnes); the other
  targets richdem (**GPL v3**). A patch that modifies GPL'd source is
  itself bound by the GPL when applied — but it is only ever applied to
  *your own* upstream clone on *your* machine for local oracle
  regeneration, and nothing is redistributed, so the obligation is
  satisfied trivially. Distributing small diffs against GPL software is
  routine and unencumbered.
- The committed test data (`test/test_cases/**`) is numeric output
  produced by the algorithms, not upstream source.

If you later set a top-level `LICENSE`, note the upstream itself is
mixed: FSM and dephier are MIT, richdem is GPL v3. That only matters if
you ever choose to *vendor* or *link* their code; this repo does
neither.
