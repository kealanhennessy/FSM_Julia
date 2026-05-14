module FillSpillMerge

# Conventions used throughout this package:
# - 1-based indexing (Julia default)
# - We index 2D arrays as M[x, y] to mirror the C++ code, even though
#   this is "row-major-like" in column-major Julia. We accept the small
#   performance cost in exchange for line-by-line translation fidelity.
# - NoData in input Float64 grids is represented as NaN.
# - In wtd outputs, ocean cells are written as 0, not NaN (matches the C++
#   fsm.exe convention and keeps the oracle comparison a plain elementwise diff).
# - See constants.jl for D8 offsets and depression sentinels;
#   fp_compare.jl for the floating-point tolerance and fp_eq/fp_le/fp_ge.

include("constants.jl")
include("fp_compare.jl")
include("bucket_fill.jl")
include("types.jl")
include("disjoint_set.jl")
include("priority_queue.jl")
include("dephier.jl")

end
