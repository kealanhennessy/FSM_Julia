module FillSpillMerge

# Conventions used throughout this package:
# - 1-based indexing (Julia default)
# - We index 2D arrays as M[x, y] to mirror the C++ code, even though
#   this is "row-major-like" in column-major Julia. We accept the small
#   performance cost in exchange for line-by-line translation fidelity.
# - NoData in Float64 grids is represented as NaN.
# - Sentinel values:
#     OCEAN     = UInt32(0)
#     NO_DEP    = typemax(UInt32)
#     NO_VALUE  = typemax(UInt32)
#     NO_FLOW   = Int8(0)
# - Floating-point tolerance for "equality" is FP_ERROR = 1e-4 (matches C++).

end
