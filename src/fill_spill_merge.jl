# Mirrors fsm/fill_spill_merge.hpp.
#
# Conventions match src/dephier.jl:
# - Matrices `topo`, `label`, `flowdirs`, `wtd` are sized `(W, H)` and
#   indexed `M[x, y]` to read like the C++.
# - Flat indices stored in `Depression` are 1-based Julia linear indices
#   (set up by dephier.jl).
# - Depression labels are 0-based (`OCEAN = 0`); look up the underlying
#   vector via `depref(deps, label)`.
# - Floating-point comparisons go through `fp_eq` / `fp_le` / `fp_ge`
#   (see fp_compare.jl). Exact `==` / `<=` are used when the C++ does so.
#
# `FP_ERROR = 1e-4` defined at the top of fill_spill_merge.hpp is dead
# code in the upstream — every `fp_*` call site resolves to the richdem
# helpers using 1e-6. We don't redefine it here.

import DataStructures

# Matches the C++ `neighbours = 8` alias; keeps loop bounds named the
# same way the reference does.
const NEIGHBOURS = 8


###############################################################################
# ResetDH
###############################################################################

"""
    reset_dh!(deps)

Reset every depression's `water_vol` to zero. Mirrors `ResetDH`.
"""
function reset_dh!(deps::DepressionHierarchy)
    for d in deps
        d.water_vol = 0.0
    end
    return
end


###############################################################################
# DepressionVolume / BackfillDepression / DetermineWaterLevel
###############################################################################

"""
    depression_volume(sill_elevation, cells_in_depression, total_elevation)

Volume held by a depression dammed at `sill_elevation`, given a cell
count and the sum of their elevations. Mirrors `DepressionVolume`.
"""
@inline depression_volume(sill_elevation, cells_in_depression::Integer, total_elevation) =
    Float64(cells_in_depression) * Float64(sill_elevation) - Float64(total_elevation)


"""
    backfill_depression!(water_level, topo, wtd, cells_affected)

Raise each affected cell's water table to the water surface elevation.
Mirrors `BackfillDepression`.
"""
function backfill_depression!(
    water_level::Float64,
    topo::AbstractMatrix,
    wtd::AbstractMatrix,
    cells_affected::AbstractVector{<:Integer},
)
    for c in cells_affected
        @assert wtd[c] == 0
        if water_level < topo[c]
            @assert fp_ge(water_level, Float64(topo[c]))
        end
        wtd[c] = water_level - Float64(topo[c])
        # Floating-point guard from the C++: never leave wtd negative
        # after a fill (this matters at saddle cells of a metadepression).
        if wtd[c] < 0
            wtd[c] = 0.0
        end
        @assert wtd[c] >= 0
    end
    return
end


"""
    determine_water_level!(wtd, cx, cy, water_vol, sill_elevation,
                           cells_to_spread_across, total_elevation) -> Float64

Lake-Level Equation. Returns the elevation of the water surface and, when
necessary, raises `wtd[cx, cy]` to absorb the residual water that won't
fit above ground. Mirrors `DetermineWaterLevel` — the C++ takes a
reference to the sill cell's wtd; we index into the array instead.
"""
function determine_water_level!(
    wtd::AbstractMatrix,
    cx::Integer,
    cy::Integer,
    water_vol::Float64,
    sill_elevation::Float64,
    cells_to_spread_across::Integer,
    total_elevation::Float64,
)
    current_dep_volume = depression_volume(sill_elevation, cells_to_spread_across, total_elevation)

    if water_vol > current_dep_volume
        if fp_eq(water_vol, current_dep_volume)
            water_vol = current_dep_volume
        end

        # The above-ground volume is exceeded but the entry condition
        # to FillDepressions guarantees this sill cell's water table can
        # absorb the rest. The assertion is the same one the C++ uses.
        sill_wtd = Float64(wtd[cx, cy])
        @assert fp_le(water_vol, current_dep_volume - sill_wtd)

        fill_amount = water_vol - current_dep_volume
        @assert fill_amount >= 0
        @assert fill_amount <= -sill_wtd
        wtd[cx, cy] = sill_wtd + fill_amount
        # Note: the C++ comment explicitly drops the `water_vol -=
        # fill_amount` update because water_vol is unused after this
        # point. We mirror that by not updating the local either.
        return sill_elevation
    elseif fp_eq(water_vol, current_dep_volume)
        return sill_elevation
    else
        # Water level below the sill: solve the Water Level Equation
        # for the level, snapping to the sill if floating-point-equal.
        nominal_level = (water_vol + total_elevation) / cells_to_spread_across
        if fp_eq(nominal_level, sill_elevation)
            return sill_elevation
        else
            return nominal_level
        end
    end
end


###############################################################################
# MoveWaterIntoPits
###############################################################################

"""
    move_water_into_pits!(topo, label, flowdirs, deps, wtd)

Route all surface water (positive `wtd`) downstream into pit cells using
the steepest-descent flow directions produced by dephier, accumulating
each leaf depression's `water_vol`. Mirrors `MoveWaterIntoPits`. All
`wtd` values are `<= 0` on exit.
"""
function move_water_into_pits!(
    topo::AbstractMatrix,
    label::AbstractMatrix,
    flowdirs::AbstractMatrix,
    deps::DepressionHierarchy,
    wtd::AbstractMatrix,
)
    W, H = size(topo)

    # Per-cell incoming-flow count. char/Int8 holds at most 8.
    dependencies = zeros(Int8, W, H)
    for y in 1:H, x in 1:W, n in 1:NEIGHBOURS
        nx = x + D8X[n]
        ny = y + D8Y[n]
        if !(1 <= nx <= W && 1 <= ny <= H)
            continue
        end
        if flowdirs[nx, ny] == D8_INVERSE[n]
            dependencies[x, y] += Int8(1)
        end
    end

    # Peaks (zero dependencies) seed the BFS. C++ iterates raw linear
    # indices in row-major order; Julia's column-major linear order over
    # an (W, H) matrix is numerically the same sequence (x varies first,
    # then y).
    q = DataStructures.Queue{Int}()
    for i in 1:length(topo)
        if dependencies[i] == 0
            push!(q, i)
        end
    end

    LI = LinearIndices(topo)
    CI = CartesianIndices(topo)

    while !isempty(q)
        c = popfirst!(q)

        # Downstream neighbour, if any. C++ uses NO_FLOW (0) as the
        # sentinel for "no downstream" both in flowdirs and as the
        # neighbour-address default.
        ndir = flowdirs[c]
        n_flat = 0   # 0 means "no downstream" (Julia 1-based; valid flat indices start at 1)
        if ndir != NO_FLOW
            ci = CI[c]
            cx, cy = ci[1], ci[2]
            nx = cx + D8X[ndir]
            ny = cy + D8Y[ndir]
            n_flat = LI[nx, ny]
            @assert n_flat >= 1
        end

        if n_flat == 0
            # Pit cell (or ocean cell — both have NO_FLOW). Push any
            # surface water into the depression's water_vol.
            if wtd[c] > 0
                depref(deps, label[c]).water_vol += Float64(wtd[c])
                wtd[c] = 0
            end
        else
            if wtd[c] > 0
                wtd[n_flat] += wtd[c]
                wtd[c] = 0
            end

            dependencies[n_flat] -= Int8(1)
            if dependencies[n_flat] == 0
                @assert dependencies[n_flat] >= 0
                push!(q, n_flat)
            end
        end
    end

    return
end


###############################################################################
# OverflowInto
###############################################################################

"""
    overflow_into!(root, stop_node, deps, jump_table, extra_water) -> dh_label_t

Recursively cascade overflowing water through the depression hierarchy.
A depression has three places to stash water — itself, its geographic
neighbour, or its parent — tried in that order. The `jump_table`
short-circuits chains of already-full depressions so the whole traversal
stays linear. Mirrors `OverflowInto`.
"""
function overflow_into!(
    root::dh_label_t,
    stop_node::dh_label_t,
    deps::DepressionHierarchy,
    jump_table::Dict{dh_label_t, dh_label_t},
    extra_water::Float64,
)
    this_dep = depref(deps, root)

    # Pick up overflow already sitting in this depression.
    if this_dep.water_vol > this_dep.dep_vol
        extra_water += this_dep.water_vol - this_dep.dep_vol
        # Exact assignment (not subtraction) for floating-point equality
        # with dep_vol, matching the C++ comment at line 527.
        this_dep.water_vol = this_dep.dep_vol
    end

    # Termination: either at the original caller's parent, or at the
    # OCEAN (which is allowed to absorb anything).
    if root == stop_node || root == OCEAN
        this_dep.water_vol += extra_water
        return root
    end

    # FIRST PLACE: in this depression.
    if this_dep.water_vol < this_dep.dep_vol
        capacity = this_dep.dep_vol - this_dep.water_vol
        if extra_water < capacity
            this_dep.water_vol = min(this_dep.water_vol + extra_water, this_dep.dep_vol)
            extra_water = 0.0
        else
            this_dep.water_vol = this_dep.dep_vol
            extra_water -= capacity
        end
    end

    if fp_eq(extra_water, 0.0)
        return root
    end

    # Jump table: shortcut through already-filled depressions.
    if haskey(jump_table, root)
        result = overflow_into!(jump_table[root], stop_node, deps, jump_table, extra_water)
        jump_table[root] = result
        return result
    end

    # SECOND PLACE: in this depression's geographic neighbour.
    if this_dep.odep != NO_VALUE
        odep = depref(deps, this_dep.odep)
        if odep.water_vol < odep.dep_vol
            result = overflow_into!(this_dep.geolink, stop_node, deps, jump_table, extra_water)
            jump_table[root] = result
            return result
        elseif odep.water_vol > odep.dep_vol
            extra_water += odep.water_vol - odep.dep_vol
            odep.water_vol = odep.dep_vol
        end
    end

    # Climb to parent. The C++ adds this and (if applicable) the
    # neighbour's water into the parent only when the parent is empty
    # and we are NOT ocean-linked — see lines 593-599.
    pdep = depref(deps, this_dep.parent)
    if pdep.water_vol == 0 && !this_dep.ocean_parent
        pdep.water_vol += this_dep.water_vol
        if this_dep.odep != NO_VALUE
            pdep.water_vol += depref(deps, this_dep.odep).water_vol
        end
    end

    # THIRD PLACE: in this depression's parent.
    result = overflow_into!(this_dep.parent, stop_node, deps, jump_table, extra_water)
    jump_table[root] = result
    return result
end


###############################################################################
# MoveWaterInDepHier
###############################################################################

"""
    move_water_in_dep_hier!(current_depression, deps, jump_table)

Post-order depth-first traversal of the depression hierarchy. After
visiting ocean-links and children, sum child volumes into the parent's
water_vol if both children are full, then trigger overflow if this
depression is overfull. Mirrors `MoveWaterInDepHier`.
"""
function move_water_in_dep_hier!(
    current_depression::dh_label_t,
    deps::DepressionHierarchy,
    jump_table::Dict{dh_label_t, dh_label_t},
)
    if current_depression == NO_VALUE
        return
    end

    this_dep = depref(deps, current_depression)

    # Visit ocean-linked depressions before children: ocean-linked can
    # flow into children but not vice-versa.
    for c in this_dep.ocean_linked
        move_water_in_dep_hier!(c, deps, jump_table)
    end

    move_water_in_dep_hier!(this_dep.lchild, deps, jump_table)
    move_water_in_dep_hier!(this_dep.rchild, deps, jump_table)

    # The ocean has no children to consolidate; bail out before we
    # mutate its water_vol.
    if current_depression == OCEAN
        return
    end

    # If both children are exactly full and the parent hasn't yet been
    # filled by an earlier OverflowInto pass, propagate their volume up.
    let lchild = this_dep.lchild, rchild = this_dep.rchild
        if lchild != NO_VALUE &&
           depref(deps, lchild).water_vol == depref(deps, lchild).dep_vol &&
           depref(deps, rchild).water_vol == depref(deps, rchild).dep_vol &&
           this_dep.water_vol == 0
            this_dep.water_vol += depref(deps, lchild).water_vol + depref(deps, rchild).water_vol
        end
    end

    if this_dep.water_vol > this_dep.dep_vol
        @assert this_dep.lchild == NO_VALUE || depref(deps, this_dep.lchild).water_vol == depref(deps, this_dep.lchild).dep_vol
        @assert this_dep.rchild == NO_VALUE || depref(deps, this_dep.rchild).water_vol == depref(deps, this_dep.rchild).dep_vol

        overflow_into!(current_depression, this_dep.parent, deps, jump_table, 0.0)

        @assert this_dep.water_vol == 0 ||
                this_dep.water_vol <= this_dep.dep_vol ||
                (this_dep.lchild == NO_VALUE && this_dep.rchild == NO_VALUE) ||
                (this_dep.lchild != NO_VALUE && this_dep.rchild != NO_VALUE &&
                 depref(deps, this_dep.lchild).water_vol < this_dep.water_vol &&
                 depref(deps, this_dep.rchild).water_vol < this_dep.water_vol)
    end

    return
end


###############################################################################
# SubtreeDepressionInfo / FindDepressionsToFill / FillDepressions
###############################################################################

"""
Tracker passed up the tree by `find_depressions_to_fill!`. Records:
- `leaf_label`: an arbitrary leaf inside the subtree (used as a pit-cell
  seed for flooding).
- `top_label`: the root of the subtree we may want to spread water over.
- `my_labels`: every depression label in the subtree (so the flood fill
  knows which cells it's allowed to enter).
"""
mutable struct SubtreeDepressionInfo
    leaf_label::dh_label_t
    top_label::dh_label_t
    my_labels::Set{dh_label_t}
end
SubtreeDepressionInfo() = SubtreeDepressionInfo(NO_VALUE, NO_VALUE, Set{dh_label_t}())


"""
    find_depressions_to_fill!(current_depression, deps, topo, label, wtd)

Depth-first descent through the hierarchy. At each node we either:
- pass our subtree info up to the parent (we're overflowing), or
- call `fill_depressions!` here and stop propagating (we have enough
  volume to contain the water).
Mirrors `FindDepressionsToFill`.
"""
function find_depressions_to_fill!(
    current_depression::dh_label_t,
    deps::DepressionHierarchy,
    topo::AbstractMatrix,
    label::AbstractMatrix,
    wtd::AbstractMatrix,
)
    if current_depression == NO_VALUE
        return SubtreeDepressionInfo()
    end

    this_dep = depref(deps, current_depression)

    # Visit ocean-linked depressions. Their water has already been
    # transferred through the metadepression tree by
    # move_water_in_dep_hier!; we discard the returned info.
    for c in this_dep.ocean_linked
        find_depressions_to_fill!(c, deps, topo, label, wtd)
    end

    if current_depression == OCEAN
        return SubtreeDepressionInfo()
    end

    left_info  = find_depressions_to_fill!(this_dep.lchild, deps, topo, label, wtd)
    right_info = find_depressions_to_fill!(this_dep.rchild, deps, topo, label, wtd)

    combined = SubtreeDepressionInfo()
    push!(combined.my_labels, current_depression)
    union!(combined.my_labels, left_info.my_labels)
    union!(combined.my_labels, right_info.my_labels)

    combined.leaf_label = left_info.leaf_label
    if combined.leaf_label == NO_VALUE
        combined.leaf_label = current_depression
    end
    combined.top_label = current_depression

    @assert this_dep.water_vol <= this_dep.dep_vol

    # Fill now if we have room, OR if our parent is ocean-linked (we
    # won't get another chance — the overflow has already left the
    # subtree), OR if our parent is dry (so we know we contain all our
    # water).
    if this_dep.water_vol < this_dep.dep_vol ||
       this_dep.ocean_parent ||
       (this_dep.water_vol == this_dep.dep_vol && depref(deps, this_dep.parent).water_vol == 0)
        @assert this_dep.water_vol <= this_dep.dep_vol

        fill_depressions!(
            depref(deps, combined.leaf_label).pit_cell,
            depref(deps, combined.top_label).out_cell,
            combined.my_labels,
            this_dep.water_vol,
            topo,
            label,
            wtd,
        )

        # Return a null marker: nothing remains to spread above us
        # (anything left has already been routed through an ocean link).
        return SubtreeDepressionInfo()
    else
        return combined
    end
end


"""
    fill_depressions!(pit_cell, out_cell, dep_labels, water_vol,
                      topo, label, wtd)

Spread `water_vol` across the metadepression containing `pit_cell`,
bounded above by `out_cell`'s elevation. Uses the
lowest-elevation-first priority queue from priority_queue.jl; on tied
elevation the most recently inserted cell pops first (matches the C++
`GridCellZk_high_pq`). Mirrors `FillDepressions`.
"""
function fill_depressions!(
    pit_cell::Integer,
    out_cell::Integer,
    dep_labels::Set{dh_label_t},
    water_vol::Float64,
    topo::AbstractMatrix,
    label::AbstractMatrix,
    wtd::AbstractMatrix,
)
    if water_vol == 0
        return
    end

    W, H = size(topo)
    LI = LinearIndices(topo)
    CI = CartesianIndices(topo)

    # Track visited cells by their flat index. Pre-size like the C++
    # reserve(2048); larger metadepressions just rehash.
    visited = Set{Int}()
    sizehint!(visited, 2048)

    flood_q = LIFOMinPriorityQueue{Tuple{Int, Int}, Float64}()

    @assert pit_cell >= 1

    pit_ci = CI[pit_cell]
    pit_x, pit_y = pit_ci[1], pit_ci[2]
    pq_push!(flood_q, Float64(topo[pit_x, pit_y]), (Int(pit_x), Int(pit_y)))
    push!(visited, Int(pit_cell))

    # Flat indices of cells whose wtd has been raised so far (used by
    # backfill_depression! once we settle on a water surface).
    cells_affected = Int[]
    total_elevation = 0.0

    out_x = CI[out_cell][1]
    out_y = CI[out_cell][2]
    out_elev = Float64(topo[out_x, out_y])

    while !isempty(flood_q)
        cx, cy = pq_pop!(flood_q)
        c_elev = Float64(topo[cx, cy])

        current_volume = depression_volume(c_elev, length(cells_affected), total_elevation)

        @assert water_vol >= 0

        if (label[cx, cy] in dep_labels) && wtd[cx, cy] > 0
            throw(ErrorException("A cell was discovered in an unfilled depression with wtd>0!"))
        end

        # Have we found a cell where the available above-ground volume
        # (plus this cell's water-table capacity) holds everything? If
        # so, settle on a water level and backfill.
        if fp_le(water_vol, current_volume - Float64(wtd[cx, cy]))
            water_level = determine_water_level!(
                wtd, cx, cy, water_vol, c_elev, length(cells_affected), total_elevation,
            )
            if fp_eq(water_level, out_elev)
                water_level = out_elev
            end

            @assert isempty(cells_affected) || fp_le(Float64(topo[cells_affected[end]]), water_level)
            @assert c_elev > water_level || fp_eq(c_elev, water_level)

            backfill_depression!(water_level, topo, wtd, cells_affected)
            return
        end

        # Not enough volume yet. Absorb this cell into the running
        # totals (skipping the outlet, which doesn't count toward volume
        # and whose wtd is the one we'd raise in DetermineWaterLevel).
        if LI[cx, cy] != out_cell
            push!(cells_affected, Int(LI[cx, cy]))
            @assert wtd[cx, cy] <= 0
            water_vol += Float64(wtd[cx, cy])   # wtd <= 0, so this can only decrease water_vol
            wtd[cx, cy] = 0
            total_elevation += c_elev
        end

        # Expand search to D8 neighbours that are part of the
        # metadepression (or are the outlet itself).
        for n in 1:NEIGHBOURS
            nx = cx + D8X[n]
            ny = cy + D8Y[n]
            if !(1 <= nx <= W && 1 <= ny <= H)
                continue
            end
            ni = LI[nx, ny]

            if !(label[ni] in dep_labels) && ni != out_cell
                continue
            end

            # Never climb past the outlet's elevation.
            if Float64(topo[nx, ny]) > out_elev
                continue
            end

            if !(ni in visited)
                pq_push!(flood_q, Float64(topo[nx, ny]), (Int(nx), Int(ny)))
                push!(visited, ni)
            end
        end

        # If the neighbour loop didn't add anything (and we'd otherwise
        # exit the outer while), fall back to the outlet so we always
        # get a final chance to settle. The C++ comments warn this can
        # infinite-loop if logic is wrong; we mirror the same fallback.
        if isempty(flood_q)
            pq_push!(flood_q, out_elev, (Int(out_x), Int(out_y)))
            push!(visited, Int(out_cell))
        end
    end

    throw(ErrorException("PQ loop exited without filling a depression!"))
end


###############################################################################
# Top-level driver
###############################################################################

"""
    fill_spill_merge!(topo, label, flowdirs, deps, wtd) -> wtd

Top-level FSM. Routes surface water into pit cells, cascades overflow
through the depression hierarchy, then spreads the resulting standing
water across the cells of each depression. Mirrors `FillSpillMerge`.

`wtd` is mutated in place; the return value is the same array for
convenience.
"""
function fill_spill_merge!(
    topo::AbstractMatrix,
    label::AbstractMatrix,
    flowdirs::AbstractMatrix,
    deps::DepressionHierarchy,
    wtd::AbstractMatrix,
)
    reset_dh!(deps)

    move_water_into_pits!(topo, label, flowdirs, deps, wtd)

    jump_table = Dict{dh_label_t, dh_label_t}()
    move_water_in_dep_hier!(dh_label_t(OCEAN), deps, jump_table)

    # Sanity checks (lines 182-187 in fill_spill_merge.hpp): every
    # non-OCEAN depression must either have no water, or hold water
    # strictly less than its volume only when its children are also
    # strictly less than its own water volume.
    for d in 1:(length(deps) - 1)
        dep = depref(deps, d)
        @assert dep.water_vol == 0 || dep.water_vol <= dep.dep_vol
        @assert dep.water_vol == 0 ||
                (dep.lchild == NO_VALUE && dep.rchild == NO_VALUE) ||
                (dep.lchild != NO_VALUE && depref(deps, dep.lchild).water_vol < dep.water_vol)
        @assert dep.water_vol == 0 ||
                (dep.lchild == NO_VALUE && dep.rchild == NO_VALUE) ||
                (dep.rchild != NO_VALUE && depref(deps, dep.rchild).water_vol < dep.water_vol)
    end

    find_depressions_to_fill!(dh_label_t(OCEAN), deps, topo, label, wtd)

    return wtd
end
