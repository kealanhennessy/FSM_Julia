# Mirrors richdem/depressions/Zhou2016.hpp — Priority-Flood (Zhou 2016
# variant). Used by Phase 4 tests as an independent reference for the
# fully-flooded hydrologic surface that FSM should produce under enough
# water for every depression to overflow.
#
# Reference: Zhou, G., Sun, Z., Fu, S., 2016. An efficient variant of the
# Priority-Flood algorithm for filling depressions in raster digital
# elevation models. Computers & Geosciences 90, 87–96.
# doi:10.1016/j.cageo.2016.02.021
#
# Indexing conventions match the rest of the package:
# - `dem` is `(W, H)` indexed `M[x, y]` to mirror C++.
# - Flat indices are 1-based Julia linear indices.
#
# The C++ uses an Int8 `labels` grid with three states: 0 (unprocessed),
# the propagated source label, and `10` for cells popped from the PQ.
# We mirror those values exactly.
#
# Priority queue entries are `(elev, flat_idx)` — mirroring the C++'s
# `std::pair<elev_t, int>` with `std::greater`. Lexicographic min: lowest
# elevation first, ties broken by lowest flat index. (The C++ tie-break
# is implementation-defined in spirit — it just falls out of pair-`<`
# semantics — but the filled-DEM output is order-invariant on the
# tie-break, so any deterministic rule produces the same result.)

import DataStructures

const _PFEntry = Tuple{Float64, Int}


# Walks "slope cells" upward from a freshly popped boundary cell.
# Mirrors `ProcessTraceQue_onepass`. Cells higher than the current focal
# cell are added to the trace queue with the focal cell's label; if the
# focal cell turns out to be a true boundary it's pushed to the priority
# queue exactly once.
function _process_trace_que_onepass!(
    dem::AbstractMatrix,
    labels::AbstractMatrix{Int8},
    trace_queue::DataStructures.Queue{Int},
    priority_queue::DataStructures.BinaryMinHeap{_PFEntry},
)
    W, H = size(dem)
    CI = CartesianIndices(dem)
    LI = LinearIndices(dem)
    while !isempty(trace_queue)
        c = popfirst!(trace_queue)
        ci = CI[c]
        cx, cy = ci[1], ci[2]
        b_in_pq = false
        for n in 1:8
            nx = cx + D8X[n]
            ny = cy + D8Y[n]
            if !(1 <= nx <= W && 1 <= ny <= H)
                continue
            end
            ni = LI[nx, ny]
            if labels[ni] != 0
                continue
            end
            if dem[c] < dem[ni]
                push!(trace_queue, ni)
                labels[ni] = labels[c]
                continue
            end

            # Decide whether c is a true boundary cell. NB: the upstream
            # Zhou2016 code uses the OUTER loop variable `n` rather than
            # the inner `nn` for the second neighbour lookup — almost
            # certainly a bug in the reference paper code. We mirror
            # exactly because: (a) the filled-DEM output is invariant
            # to which sub-check classifies a cell as a boundary
            # (it just changes PQ traffic), and (b) deviating would
            # mean we no longer match the C++ at all on traversal
            # statistics. If the test ever asserts on PQ size, we'd
            # need to revisit.
            if !b_in_pq
                is_boundary = true
                for _nn in 1:8
                    nnx = nx + D8X[n]
                    nny = ny + D8Y[n]
                    if !(1 <= nnx <= W && 1 <= nny <= H)
                        continue
                    end
                    nni = LI[nnx, nny]
                    if labels[nni] != 0 && dem[nni] < dem[ni]
                        is_boundary = false
                        break
                    end
                end
                if is_boundary
                    push!(priority_queue, (Float64(dem[c]), c))
                    b_in_pq = true
                end
            end
        end
    end
    return
end


# Drain a depression by raising its cells to the spill elevation. Mirrors
# `ProcessPit_onepass`. Neighbours above the dam go to the trace queue;
# at-or-below cells are raised to `c_elev` and re-queued.
function _process_pit_onepass!(
    c_elev::Float64,
    dem::AbstractMatrix,
    labels::AbstractMatrix{Int8},
    depression_que::DataStructures.Queue{Int},
    trace_queue::DataStructures.Queue{Int},
    priority_queue::DataStructures.BinaryMinHeap{_PFEntry},
)
    W, H = size(dem)
    CI = CartesianIndices(dem)
    LI = LinearIndices(dem)
    while !isempty(depression_que)
        c = popfirst!(depression_que)
        ci = CI[c]
        cx, cy = ci[1], ci[2]
        for n in 1:8
            nx = cx + D8X[n]
            ny = cy + D8Y[n]
            if !(1 <= nx <= W && 1 <= ny <= H)
                continue
            end
            ni = LI[nx, ny]
            if labels[ni] != 0
                continue
            end
            labels[ni] = labels[c]
            if Float64(dem[ni]) > c_elev
                push!(trace_queue, ni)
            else
                dem[ni] = c_elev
                push!(depression_que, ni)
            end
        end
    end
    return
end


"""
    priority_flood_zhou2016!(dem) -> dem

Fill all pits and remove all digital dams from `dem`, in place. Mirrors
`PriorityFlood_Zhou2016`. Boundary cells (the outer ring) are treated as
the spill points; depressions are raised to the level of their lowest
boundary path. Returns the mutated `dem`.
"""
function priority_flood_zhou2016!(dem::AbstractMatrix)
    W, H = size(dem)

    trace_queue    = DataStructures.Queue{Int}()
    depression_que = DataStructures.Queue{Int}()

    labels = zeros(Int8, W, H)
    priority_queue = DataStructures.BinaryMinHeap{_PFEntry}()

    LI = LinearIndices(dem)
    @inline place_cell!(x, y) = push!(priority_queue, (Float64(dem[x, y]), LI[x, y]))

    # Seed the outer ring. Top + bottom rows, then left + right columns
    # excluding corners (matches the C++ four-loop layout).
    for x in 1:W; place_cell!(x, 1); end
    for x in 1:W; place_cell!(x, H); end
    for y in 2:(H - 1); place_cell!(1, y); end
    for y in 2:(H - 1); place_cell!(W, y); end

    CI = CartesianIndices(dem)
    while !isempty(priority_queue)
        elev_seed, c = pop!(priority_queue)
        labels[c] = Int8(10)
        ci = CI[c]
        cx, cy = ci[1], ci[2]

        for n in 1:8
            nx = cx + D8X[n]
            ny = cy + D8Y[n]
            if !(1 <= nx <= W && 1 <= ny <= H)
                continue
            end
            ni = LI[nx, ny]
            if labels[ni] != 0
                continue
            end
            labels[ni] = labels[c]
            if Float64(dem[ni]) <= elev_seed
                dem[ni] = elev_seed
                push!(depression_que, ni)
                _process_pit_onepass!(elev_seed, dem, labels, depression_que, trace_queue, priority_queue)
            else
                push!(trace_queue, ni)
            end
            _process_trace_que_onepass!(dem, labels, trace_queue, priority_queue)
        end
    end

    return dem
end
