# Case 4: Nested metadepression

10×10 grid. Three-level depression hierarchy:

- Outer plateau (elev 9) surrounded by NoData ocean.
- A "big bowl" inside the plateau: rim at elev 5, floor at elev 3.
- Two small pits inside the big bowl's floor:
  - Pit A (rows 4–5, col 3): elev 1
  - Pit B (rows 4–5, col 6): elev 2

Pits A and B are separated by elev-3 floor cells, so they are distinct
depressions at the lowest level of the hierarchy. Both spill at elev 3
into the surrounding big bowl. The big bowl spills at elev 5 over its
rim into the plateau (and from there to the ocean).

## What this exercises

The full depression hierarchy traversal: water rises in A and B until
they spill at elev 3 into the big bowl, the bowl then accumulates water
across both pits plus its own floor cells, eventually all three depression
levels are at a common water level if filling continues past elev 3.

## Input

- `ocean_level = 0`
- `--swl 1.0`

## Expected output (`expected-wtd.tif`)

Common water level in the merged structure: **4.875** (just below the
big bowl rim at elev 5).

- Pit A cells (elev 1): wtd = 3.875
- Pit B cells (elev 2): wtd = 2.875
- Big bowl floor (elev 3): wtd = 1.875
- Bowl rim (elev 5), plateau (elev 9), ocean: wtd = 0

Volume check: 2×3.875 + 2×2.875 + 12×1.875 = 36 units in the merged
structure (out of 64 total input; the remaining 28 drained to ocean from
the plateau and outer rim).
