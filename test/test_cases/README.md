# FSM test cases for the Julia port

Input/output pairs from the reference C++ Fill-Spill-Merge binary, for
verifying a Julia port of the algorithm against the upstream
implementation.

## Layout

Each `case_NN_*/` directory contains:

- `input.asc` — Arc/Info ASCII grid of the DEM (human-readable source of truth).
- `input.tif` — same DEM converted to a Float64 GeoTIFF (the format FSM consumes).
- `expected-wtd.tif` — water table depth produced by the C++ binary
  (water depth above the topographic surface, in input units).
- `expected-hydrologic-surface-height.tif` — `wtd + topo`, i.e. the
  height of the water surface where it exists, or topography elsewhere.
- `run.sh` — reproduces the test from `input.asc`.
- `README.md` — description of what the case exercises and the expected
  output structure.

Case 5 also has `generate.py` for re-deriving the 100×100 input.

## Algorithm contract

The C++ binary `build/fsm.exe` takes:

```
fsm.exe <topography.tif> <output_prefix> <ocean_level> [--swl <X> | --swf <file>]
```

- **`ocean_level`**: a scalar elevation. FSM bucket-fills from the grid
  edges inward, marking any edge-connected cell with elev ≤ `ocean_level`
  as OCEAN. **Additionally, any NoData cell becomes OCEAN regardless of
  position.** FSM errors out if no OCEAN cells are found, so every test
  case here has at least one NoData cell.
- **`--swl <X>`**: constant surface water level — every non-ocean cell
  starts with X units of water on it. (`--swf <file>` is the alternative
  for non-uniform water input; not used in these test cases.)
- **Output**: two GeoTIFFs at `<output_prefix>-wtd.tif` (water table
  depth) and `<output_prefix>-hydrologic-surface-height.tif` (wtd + topo).
- **Ocean cells**: have `wtd = 0` in the output.

## Cases

| # | Size    | Purpose                                                              |
|---|---------|----------------------------------------------------------------------|
| 1 | 5×5     | Trivial trough, one NoData ocean cell. Sanity check.                 |
| 2 | 10×10   | Single bowl-shaped depression, NoData border ocean. No overflow.     |
| 3 | 10×10   | Two adjacent depressions sharing a sill. Forces overflow + merge.    |
| 4 | 10×10   | Three-level hierarchy: two pits inside a big bowl inside plateau.    |
| 5 | 100×100 | Synthetic 1/f-noise terrain. Many depressions; full-scale exercise.  |

## Regenerating

To regenerate every output from scratch (e.g. after rebuilding `fsm.exe`):

```
./run_all.sh
```

Each individual case can also be re-run via its own `./run.sh`. The path
to the FSM binary can be overridden via the `FSM_EXE` environment variable.

## Using these from Julia

For each case, read `input.tif`, run the Julia port with the same
`ocean_level` and `--swl` values listed in that case's `README.md` and
`run.sh`, then compare the Julia output against `expected-wtd.tif`
cell-by-cell. Floating-point tolerance of ~1e-9 should be safe given the
algorithm is fundamentally arithmetic on volumes; tighter is fine for
cases 1–4 where outputs are exact rationals.
