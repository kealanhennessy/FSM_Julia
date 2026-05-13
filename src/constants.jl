# Constants ported from the C++ reference (originals in
# richdem/common/constants.hpp and dephier/dephier.hpp).

# D8 neighbour offsets. The C++ arrays have 9 elements with index 0 unused;
# here we drop that placeholder so direction codes 1..8 index the tuples
# directly.
# Layout (direction code at each position relative to centre `.`):
#
#     2 3 4       NW N  NE
#     1 . 5       W  .  E
#     8 7 6       SW S  SE
#
# x increases east, y increases south.
const D8X = (-1, -1,  0,  1, 1, 1, 0, -1)
const D8Y = ( 0, -1, -1, -1, 0, 1, 1,  1)
# D8_INVERSE[n] = direction code pointing from neighbour n back to the centre.
const D8_INVERSE = (5, 6, 7, 8, 1, 2, 3, 4)

const NO_FLOW = Int8(0)

# Label-grid sentinels. OCEAN tags ocean cells in the label Matrix{UInt32};
# NO_DEP tags cells before depression assignment. NO_VALUE and NO_PARENT are
# the same typemax(UInt32) value, kept as distinct names to match the C++
# usage at the relevant call sites (cell indices, parent links).
const OCEAN     = UInt32(0)
const NO_DEP    = typemax(UInt32)
const NO_VALUE  = typemax(UInt32)
const NO_PARENT = typemax(UInt32)
