#!/usr/bin/env python3
"""Generate a deterministic 100x100 pseudo-Perlin DEM as an Arc/Info ASCII grid.

The DEM is built from 1/f spectral noise (a common "pink noise" trick that
produces natural-looking terrain), then a 1-cell NoData border is added so
FSM sees an ocean ring around the entire grid.
"""
import numpy as np

N = 100
SEED = 42
NODATA = -9999

rng = np.random.default_rng(SEED)

# 1/f-style spectral noise: random Fourier coefficients weighted by 1/k^alpha.
freq = np.fft.fftfreq(N)
fx, fy = np.meshgrid(freq, freq, indexing="xy")
k = np.sqrt(fx**2 + fy**2)
k[0, 0] = 1.0  # avoid div-by-zero at DC

amplitude = 1.0 / (k ** 1.8)
amplitude[0, 0] = 0.0  # zero out DC

phase = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
elev = np.real(np.fft.ifft2(phase * amplitude))

# Rescale to 0..50 metre-ish elevations.
elev -= elev.min()
elev /= elev.max()
elev *= 50.0

# Punch out a NoData border so FSM has ocean cells.
elev[0, :] = NODATA
elev[-1, :] = NODATA
elev[:, 0] = NODATA
elev[:, -1] = NODATA

header = (
    f"ncols        {N}\n"
    f"nrows        {N}\n"
    f"xllcorner    0\n"
    f"yllcorner    0\n"
    f"cellsize     1\n"
    f"NODATA_value {NODATA}\n"
)

with open("input.asc", "w") as f:
    f.write(header)
    for row in elev:
        f.write(" ".join(f"{v:.6f}" if v != NODATA else str(NODATA) for v in row))
        f.write("\n")

print(f"Wrote input.asc (seed={SEED}, range {elev[elev != NODATA].min():.3f} .. {elev.max():.3f})")
