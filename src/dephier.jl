# Mirrors GetDepressionHierarchy + CalculateMarginalVolumes +
# CalculateTotalVolumes from
# submodules/dephier/include/dephier/dephier.hpp.
#
# Single-threaded (the C++ uses OpenMP for the seed-finding and marginal-
# volume passes; we just do the work sequentially since dephier has only
# ever been run on small-to-medium grids in this project).
#
# Indexing conventions:
# - Grids `dem`, `label`, `flowdirs` are size `(W, H)` and indexed
#   `M[x, y]` to mirror C++.
# - Flat indices stored in `Depression`/`Outlet` are 1-based Julia
#   linear indices. The dump-file loader converts the 0-based C++
#   indices on read; nothing inside the algorithm needs to know about
#   the 0/1 distinction.
# - Depression labels are 0-based (OCEAN = 0). To index into the
#   `depressions::Vector{Depression}`, add 1: `depressions[Int(label)+1]`.

# Helper: the Julia `depressions` vector is 1-based but labels are
# 0-based; this hides the +1.
@inline depref(deps::DepressionHierarchy, label::Integer) = deps[Int(label) + 1]


"""
    get_depression_hierarchy(dem, label, flowdirs) -> DepressionHierarchy{T}

Compute the depression hierarchy for `dem`. `label` must be pre-populated
with `OCEAN` for ocean cells and `NO_DEP` everywhere else. Both `label`
and `flowdirs` are mutated in place: labels are filled in for every
non-ocean cell, and flowdirs are set for every cell that has a downhill
neighbour (pit cells keep `NO_FLOW`).

The returned vector contains the leaf depressions (one per pit cell)
followed by the meta-depressions formed by merging children when their
shared outlet is processed. `depressions[1]` is always the OCEAN
depression (label 0).
"""
function get_depression_hierarchy(
    dem::AbstractMatrix{T},
    label::AbstractMatrix{dh_label_t},
    flowdirs::AbstractMatrix{Int8},
) where {T<:AbstractFloat}
    W, H = size(dem)
    size(label)    == (W, H) || throw(DimensionMismatch("label size must match dem"))
    size(flowdirs) == (W, H) || throw(DimensionMismatch("flowdirs size must match dem"))

    depressions = DepressionHierarchy{T}()

    outlet_database = Dict{OutletLink, Outlet{T}}()

    LI = LinearIndices(label)
    CI = CartesianIndices(label)

    # ---- Pass 1: collect ocean seeds (ocean cells with at least one
    # non-ocean neighbour) and validate input labels.
    ocean_seeds = flat_c_idx[]
    sizehint!(ocean_seeds, length(label) ÷ 40)
    ocean_cells = 0

    for y in 1:H, x in 1:W
        if label[x, y] != OCEAN
            if label[x, y] != NO_DEP
                throw(ArgumentError(
                    "Label array given to get_depression_hierarchy must contain only NO_DEP and OCEAN labels!",
                ))
            end
            continue
        end
        has_non_ocean = false
        for n in 1:8
            nx = x + D8X[n]
            ny = y + D8Y[n]
            if 1 <= nx <= W && 1 <= ny <= H && label[nx, ny] != OCEAN
                has_non_ocean = true
                break
            end
        end
        if has_non_ocean
            push!(ocean_seeds, flat_c_idx(LI[x, y]))
            ocean_cells += 1
        end
    end

    if ocean_cells == 0
        throw(ArgumentError("No OCEAN cells found, could not make a DepressionHierarchy!"))
    end

    # ---- Initialize the OCEAN depression at index 0 (Julia slot 1).
    oceandep = Depression{T}()
    oceandep.pit_elev = T(-Inf)
    oceandep.pit_cell = NO_VALUE
    oceandep.dep_label = dh_label_t(0)
    push!(depressions, oceandep)

    # ---- Pass 2: collect pit cells (no lower neighbour).
    land_seeds = flat_c_idx[]
    sizehint!(land_seeds, length(label) ÷ 40)
    pit_cell_count = 0

    for y in 1:H, x in 1:W
        if label[x, y] == OCEAN
            continue
        end
        my_elev = dem[x, y]
        has_lower = false
        for n in 1:8
            nx = x + D8X[n]
            ny = y + D8Y[n]
            if !(1 <= nx <= W && 1 <= ny <= H)
                continue
            end
            if dem[nx, ny] < my_elev
                has_lower = true
                break
            end
        end
        if !has_lower
            push!(land_seeds, flat_c_idx(LI[x, y]))
            pit_cell_count += 1
        end
    end

    # Mirror C++: sort seeds for determinism (they're collected in y-major
    # order so already sorted, but keep the sort for parity with the
    # parallel C++ path that reduces from threads).
    sort!(ocean_seeds)
    sort!(land_seeds)

    # ---- Pass 3: priority-queue traversal. Pop lowest elev first; for
    # equal elevations, most recently inserted comes first (LIFO).
    pq = LIFOMinPriorityQueue{flat_c_idx, T}()
    for ci in ocean_seeds
        pq_push!(pq, dem[ci], ci)
    end
    empty!(ocean_seeds); sizehint!(ocean_seeds, 0)
    for ci in land_seeds
        pq_push!(pq, dem[ci], ci)
    end
    empty!(land_seeds); sizehint!(land_seeds, 0)

    sizehint!(outlet_database, 3 * (pit_cell_count + 1))

    while !isempty(pq)
        ci    = pq_pop!(pq)
        celev = dem[ci]
        clabel = label[ci]
        cxy = CI[ci]
        cx, cy = cxy[1], cxy[2]

        if clabel == OCEAN
            # Already in queue; nothing to do.
        elseif clabel == NO_DEP
            # First time we've reached this pit cell; create a new
            # depression for it.
            clabel = dh_label_t(length(depressions))
            newdep = Depression{T}()
            newdep.pit_cell  = flat_c_idx(ci)
            newdep.pit_elev  = celev
            newdep.dep_label = clabel
            push!(depressions, newdep)
            label[ci] = clabel
        else
            # Cell on the frontier of an existing depression — fall
            # through and inspect its neighbours.
        end

        for n in 1:8
            nx = cx + D8X[n]
            ny = cy + D8Y[n]
            if !(1 <= nx <= W && 1 <= ny <= H)
                continue
            end
            ni = LI[nx, ny]
            nlabel = label[ni]

            if nlabel == NO_DEP
                label[ni] = clabel
                pq_push!(pq, dem[ni], flat_c_idx(ni))
                flowdirs[nx, ny] = D8_INVERSE[n]
            elseif nlabel == clabel
                # Same depression; nothing to do (common case for flats).
            else
                # Found neighbouring depression. The outlet between the
                # two is whichever of the focal/neighbour cells is higher.
                out_cell = flat_c_idx(ci)
                out_elev = celev
                if dem[ni] > out_elev
                    out_cell = flat_c_idx(ni)
                    out_elev = dem[ni]
                end

                # Note: clabel/nlabel order is intentionally NOT swapped
                # here — we mirror the C++ which keeps both (a,b) and
                # (b,a) keys; both end up with identical out_elev because
                # the candidate outlet cell is symmetric.
                olink = (clabel, nlabel)
                existing = get(outlet_database, olink, nothing)
                if existing === nothing
                    outlet_database[olink] = Outlet{T}(clabel, nlabel, out_cell, out_elev)
                elseif existing.out_elev > out_elev
                    existing.out_cell = out_cell
                    existing.out_elev = out_elev
                end
            end
        end
    end

    # ---- Pass 4: copy outlets out of the hash table into a vector and
    # sort by elevation, then build the hierarchy.
    outlets = Vector{Outlet{T}}(undef, length(outlet_database))
    let i = 1
        for o in values(outlet_database)
            outlets[i] = o
            i += 1
        end
    end
    empty!(outlet_database)

    # Sort by (out_elev, depa, depb, out_cell) for full determinism. The
    # dephier algorithm itself is insensitive to which equally-low outlet
    # is picked (any choice yields a hierarchy that captures the same
    # physical depressions), but a deterministic tie-break makes the
    # Julia↔C++ oracle comparison bit-exact. The patched dephier.hpp
    # uses the same comparator.
    sort!(outlets, by = o -> (o.out_elev, o.depa, o.depb, o.out_cell))

    djset = DisjointDenseIntSet(length(depressions))

    for outlet in outlets
        depa_set = find_set!(djset, outlet.depa)
        depb_set = find_set!(djset, outlet.depb)

        if depa_set == depb_set
            continue
        end

        if depa_set == OCEAN || depb_set == OCEAN
            # Ensure `depb` is the ocean side.
            if depa_set == OCEAN
                outlet.depa, outlet.depb = outlet.depb, outlet.depa
                depa_set, depb_set = depb_set, depa_set
            end

            dep = depref(depressions, depa_set)
            @assert dep.out_cell == NO_VALUE
            @assert dep.odep     == NO_VALUE

            dep.parent       = outlet.depb
            dep.out_elev     = outlet.out_elev
            dep.out_cell     = outlet.out_cell
            dep.odep         = NO_VALUE
            dep.ocean_parent = true
            dep.geolink      = outlet.depb
            push!(depref(depressions, outlet.depb).ocean_linked, depa_set)
            merge_a_into_b!(djset, depa_set, OCEAN)
        else
            depa = depref(depressions, depa_set)
            depb = depref(depressions, depb_set)
            @assert depa.odep == NO_VALUE
            @assert depb.odep == NO_VALUE

            newlabel = dh_label_t(length(depressions))
            depa.parent   = newlabel
            depb.parent   = newlabel
            depa.out_cell = outlet.out_cell
            depb.out_cell = outlet.out_cell
            depa.out_elev = outlet.out_elev
            depb.out_elev = outlet.out_elev
            depa.odep     = depb_set
            depb.odep     = depa_set
            depa.geolink  = outlet.depb
            depb.geolink  = outlet.depa

            depa_pitcell_temp = depa.pit_cell

            newdep = Depression{T}()
            newdep.lchild    = depa_set
            newdep.rchild    = depb_set
            newdep.dep_label = newlabel
            newdep.pit_cell  = depa_pitcell_temp
            push!(depressions, newdep)

            merge_a_into_b!(djset, depa_set, newlabel)
            merge_a_into_b!(djset, depb_set, newlabel)
        end
    end

    calculate_marginal_volumes!(depressions, dem, label)
    calculate_total_volumes!(depressions)

    return depressions
end


"""
    calculate_marginal_volumes!(deps, dem, label)

Walks each cell upward through its chain of meta-parents until it finds
the lowest meta-depression whose `out_elev` exceeds the cell's elevation,
then accumulates the cell's count and elevation into that depression's
marginal totals (`cell_count`, `total_elevation`). Cells that walk all
the way to the ocean contribute nothing.
"""
function calculate_marginal_volumes!(
    deps::DepressionHierarchy{T},
    dem::AbstractMatrix{T},
    label::AbstractMatrix{dh_label_t},
) where {T<:AbstractFloat}
    for i in eachindex(dem)
        my_elev = dem[i]
        clabel  = label[i]

        while clabel != OCEAN && my_elev > depref(deps, clabel).out_elev
            clabel = depref(deps, clabel).parent
        end

        if clabel == OCEAN
            continue
        end

        d = depref(deps, clabel)
        d.cell_count       += UInt32(1)
        d.total_elevation  += dem[i]
    end
    return
end


"""
    calculate_total_volumes!(deps)

Propagates per-leaf cell counts and elevations up the binary tree, then
computes each depression's `dep_vol` from the accumulated totals. Relies
on the invariant that child labels are always less than their parent's
label (established during hierarchy construction).
"""
function calculate_total_volumes!(deps::DepressionHierarchy{T}) where {T<:AbstractFloat}
    for d in 0:length(deps)-1
        dep = depref(deps, d)
        if dep.lchild != NO_VALUE
            @assert dep.rchild != NO_VALUE
            @assert dep.lchild < d
            @assert dep.rchild < d
            lc = depref(deps, dep.lchild)
            rc = depref(deps, dep.rchild)
            dep.cell_count      += lc.cell_count + rc.cell_count
            dep.total_elevation += lc.total_elevation + rc.total_elevation
        end
        # NB: `out_elev` is `T(Inf)` for the OCEAN and any unmerged
        # top-level depression that never found an outlet. The product
        # `0 * Inf` is `NaN`; the C++ produces the same value.
        dep.dep_vol = Float64(dep.cell_count) * Float64(dep.out_elev) - dep.total_elevation
        @assert dep.lchild == NO_VALUE || fp_le(
            depref(deps, dep.lchild).dep_vol + depref(deps, dep.rchild).dep_vol,
            dep.dep_vol,
        )
    end
    return
end
