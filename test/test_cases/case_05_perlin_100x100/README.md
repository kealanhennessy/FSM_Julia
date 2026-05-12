# Case 5: 100×100 pseudo-Perlin DEM

Synthetic 100×100 "natural-looking" terrain generated deterministically
from a 1/f spectral noise (similar in character to Perlin noise — many
local minima, smooth gradients). Edges are NoData → ocean.

I went with synthetic terrain instead of real elevation data because the
goal is reproducibility for the Julia port: anyone can regenerate the
exact same input from a fixed seed, no external downloads required. If
you want a real SRTM tile instead, drop the GeoTIFF in as `input.tif`
and re-run `run.sh` (skipping the Python generation step).

## What this exercises

- A non-trivial number of depressions (~hundreds of pits and
  metadepressions) with a real depression hierarchy.
- A 100×100 grid — large enough that subtle correctness bugs in the
  Julia port should show up as numerical divergence somewhere.
- The full pipeline at realistic scale, not just hand-crafted edge cases.

## Input

- `generate.py` — produces `input.asc` deterministically (`seed = 42`).
- Elevation range: 0 .. 50, with NoData border one cell thick.
- `ocean_level = 0`, `--swl 0.5`.

## Expected output (`expected-wtd.tif`)

- ~979 of 9604 non-ocean cells receive water (>0.01 units deep).
- Max wtd ≈ 4.9 units.
- Total water volume in depressions ≈ 1059 units (out of 9604 × 0.5 ≈ 4802
  total input; the rest drained to ocean).

The Julia port should reproduce these aggregate stats and match
`expected-wtd.tif` cell-for-cell within floating-point tolerance.

## Reproduce

```
./run.sh
```

(Requires Python 3 with numpy.)
