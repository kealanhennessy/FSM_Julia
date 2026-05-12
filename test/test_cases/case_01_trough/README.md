# Case 1: Trivial trough

5×5 grid. Middle row at elevation 1, all other rows at elevation 10. One
NoData cell at (0,0) acts as the single OCEAN cell (FSM requires at least one).

## What this exercises

- Basic depression filling — water on high cells drains down into the trough.
- A single isolated ocean cell that receives some drainage from the
  immediately adjacent high cells but is otherwise disconnected from the
  trough.

## Input

```
NoData  10  10  10  10
   10   10  10  10  10
    1    1   1   1   1
   10   10  10  10  10
   10   10  10  10  10
```

- `ocean_level = -100` (only the NoData cell becomes OCEAN)
- `--swl 1.0` (every non-ocean cell starts with 1.0 units of water)

## Expected output (`expected-wtd.tif`)

Trough cells (row 2) hold ≈4.2 units each; all other non-ocean cells are 0.
Ocean cell is 0. The "missing" ~3 units (vs. naive 24/5 = 4.8) drained
to the corner ocean via the adjacent plateau cells.

## Reproduce

```
./run.sh
```
