# Case 3: Two adjacent depressions sharing a sill

10×10 grid. Edge NoData → ocean. Two depressions side-by-side at
different floor elevations, separated by a one-column sill:

- Depression A (cols 2–3): elev 1 (lower)
- Sill (col 4): elev 3
- Depression B (cols 5–6): elev 2 (higher)

Outer plateau at elev 9 surrounds both depressions.

## What this exercises

The merge step of FSM. With enough water, the lower depression A fills
to the sill, overflows into B, and the two become a single combined
depression filled to a common water level. This test verifies that
overflow/merge accounting is correct.

(With smaller `--swl` neither overflows — see commit history if you want
to also test the "separate depressions" regime; here we drive enough
water to force the merge.)

## Input

- `ocean_level = 0`
- `--swl 1.5`

## Expected output (`expected-wtd.tif`)

A, sill, and B are all submerged at common water level **4.5**:

- A cells (elev 1): wtd = 3.5
- Sill cells (elev 3): wtd = 1.5
- B cells (elev 2): wtd = 2.5
- Plateau and ocean: wtd = 0

Volume check: 8×3.5 + 4×1.5 + 8×2.5 = 54 units in the merged depression
(out of 96 total input; the remaining 42 drained to ocean from the
outer plateau).
