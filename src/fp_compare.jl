# Floating-point comparison helpers. Mirrors richdem/common/math.hpp.
#
# fsm/fill_spill_merge.hpp also defines `FP_ERROR = 1e-4`, but it's dead in
# the C++ source — every fp_eq/fp_le/fp_ge call site resolves to the richdem
# implementation below, which uses 1e-6.
const FP_COMPARISON_ERROR = 1e-6

fp_le(a::Real, b::Real) = a < b || abs(a - b) < FP_COMPARISON_ERROR
fp_ge(a::Real, b::Real) = a > b || abs(a - b) < FP_COMPARISON_ERROR
fp_eq(a::Real, b::Real) = abs(a - b) < FP_COMPARISON_ERROR
