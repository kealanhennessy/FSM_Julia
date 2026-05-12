# Case 2: Single depression with ocean

10×10 grid. NoData border on all four sides forms the ocean. Interior is
a plateau at elevation 5 with a single bowl in the middle: a 2×2 center
at elevation 1, surrounded by a one-cell rim at elevation 3.

## What this exercises

- Ocean detected from edge NoData cells (the canonical FSM setup).
- A single depression whose sill (elev 3) is below the surrounding plateau
  (elev 5), so plateau cells closer to the bowl drain inward, and plateau
  cells closer to the edge drain to the ocean.
- The depression does not overflow into the ocean — water pools at a
  level below the rim of the surrounding plateau.

## Input

- Edges (row 0, row 9, col 0, col 9): NoData → OCEAN
- Plateau: elev 5
- Bowl rim (12 cells at rows 3–6, cols 3–6 minus center): elev 3
- Bowl center (4 cells at rows 4–5, cols 4–5): elev 1
- `ocean_level = 0`, `--swl 1.0`

## Expected output (`expected-wtd.tif`)

Water level in the bowl rises to ≈4.75:
- Center cells (elev 1): wtd = 3.75
- Rim cells (elev 3): wtd = 1.75
- Plateau: wtd = 0 (drained to ocean)
- Ocean: wtd = 0

Volume check: 4×3.75 + 12×1.75 = 36 units, which equals 16 bowl cells' own
water (16) plus 20 units of plateau drainage. The other 28 plateau cells
drained to the ocean.
