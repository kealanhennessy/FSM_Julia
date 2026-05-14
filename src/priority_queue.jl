# Min-heap with LIFO behaviour for equal keys, mirroring the C++
# `radix_heap::pair_radix_heap` used by dephier.
#
# Why LIFO matters: the dephier algorithm relies on a flood-fill where
# all cells of equal elevation in a single flat area should end up in
# one depression. The C++ pulls equal-elev items in reverse-insertion
# order so the wavefront keeps consuming the *current* flat before
# revisiting the original seed cell — see the comment block at
# dephier.hpp:413-416.
#
# Why a BinaryHeap rather than DataStructures.PriorityQueue: the
# algorithm legitimately pushes the same cell twice when a pit-seed
# cell is also reached by the watershed wavefront before being popped.
# PriorityQueue is dict-backed and rejects duplicate keys, so we use
# a plain heap of `(key, -counter, value)` tuples; tuple comparison
# gives both the elev ordering and the LIFO tie-break for free.

import DataStructures

mutable struct LIFOMinPriorityQueue{V, K}
    heap::DataStructures.BinaryMinHeap{Tuple{K, Int, V}}
    counter::Int
end

function LIFOMinPriorityQueue{V, K}() where {V, K}
    return LIFOMinPriorityQueue{V, K}(
        DataStructures.BinaryMinHeap{Tuple{K, Int, V}}(),
        0,
    )
end

Base.length(q::LIFOMinPriorityQueue)  = length(q.heap)
Base.isempty(q::LIFOMinPriorityQueue) = isempty(q.heap)

function pq_push!(q::LIFOMinPriorityQueue{V, K}, key::K, value::V) where {V, K}
    q.counter += 1
    push!(q.heap, (key, -q.counter, value))
    return q
end

# Returns the value (e.g. a flat cell index). Caller can re-derive the
# key from whatever array the value indexes into (e.g. `dem[ci]`).
function pq_pop!(q::LIFOMinPriorityQueue)
    return pop!(q.heap)[3]
end
