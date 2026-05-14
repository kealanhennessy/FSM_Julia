# Helpers for the Phase 2 dephier oracle tests.
#
# Parses the `expected-dh.txt` files produced by the C++ dh_dump tool
# (see Barnes2020-FillSpillMerge worktree `tools/dh_dump.cpp`). Format:
#
#   WIDTH <W>
#   HEIGHT <H>
#   NDEP <N>
#   LABEL
#   <H rows of W ints>
#   FLOWDIRS
#   <H rows of W ints>
#   DEPRESSIONS
#   # column header
#   <N rows>
#
# Flat indices in the C++ dump are 0-based; the parser converts to
# 1-based on read so the resulting fields can be compared directly
# against Julia `Depression` structs.

import FillSpillMerge:
    OCEAN, NO_DEP, NO_FLOW, NO_VALUE,
    bucket_fill_from_edges!,
    Depression, dh_label_t, flat_c_idx

# Parser-friendly variant of `parse(Float64, s)` that handles the
# lowercase `inf`/`nan` C++ ostream produces across Julia versions.
function _parse_float(s::AbstractString)
    sl = lowercase(s)
    if sl == "nan" || sl == "-nan"
        return NaN
    elseif sl == "inf" || sl == "+inf"
        return Inf
    elseif sl == "-inf"
        return -Inf
    else
        return parse(Float64, s)
    end
end

# 0-based C++ flat-index -> 1-based Julia linear index, preserving the
# NO_VALUE sentinel.
_to_julia_flat(i::UInt32) = i == NO_VALUE ? NO_VALUE : (i + UInt32(1))

struct OracleDH
    W::Int
    H::Int
    label::Matrix{dh_label_t}
    flowdirs::Matrix{Int8}
    deps::Vector{Depression{Float64}}
end

function read_dh_dump(path::AbstractString)
    lines = readlines(path)
    W = parse(Int, split(lines[1])[2])
    H = parse(Int, split(lines[2])[2])
    N = parse(Int, split(lines[3])[2])

    @assert lines[4] == "LABEL"
    label = Matrix{dh_label_t}(undef, W, H)
    for y in 1:H
        row = parse.(dh_label_t, split(lines[4 + y]))
        @assert length(row) == W
        for x in 1:W
            label[x, y] = row[x]
        end
    end

    fd_hdr_idx = 4 + H + 1
    @assert lines[fd_hdr_idx] == "FLOWDIRS"
    flowdirs = Matrix{Int8}(undef, W, H)
    for y in 1:H
        row = parse.(Int, split(lines[fd_hdr_idx + y]))
        @assert length(row) == W
        for x in 1:W
            flowdirs[x, y] = row[x]
        end
    end

    deps_hdr_idx = fd_hdr_idx + H + 1
    @assert lines[deps_hdr_idx] == "DEPRESSIONS"
    @assert startswith(lines[deps_hdr_idx + 1], "#")

    deps = Vector{Depression{Float64}}(undef, N)
    for i in 1:N
        cols = split(lines[deps_hdr_idx + 1 + i])
        # Columns: idx pit_cell out_cell parent odep geolink pit_elev
        # out_elev lchild rchild ocean_parent dep_label cell_count
        # total_elevation dep_vol water_vol ocean_linked
        d = Depression{Float64}()
        d.pit_cell        = _to_julia_flat(parse(UInt32, cols[2]))
        d.out_cell        = _to_julia_flat(parse(UInt32, cols[3]))
        d.parent          = parse(UInt32, cols[4])
        d.odep            = parse(UInt32, cols[5])
        d.geolink         = parse(UInt32, cols[6])
        d.pit_elev        = _parse_float(cols[7])
        d.out_elev        = _parse_float(cols[8])
        d.lchild          = parse(UInt32, cols[9])
        d.rchild          = parse(UInt32, cols[10])
        d.ocean_parent    = parse(Int, cols[11]) != 0
        d.dep_label       = parse(UInt32, cols[12])
        d.cell_count      = parse(UInt32, cols[13])
        d.total_elevation = _parse_float(cols[14])
        d.dep_vol         = _parse_float(cols[15])
        d.water_vol       = _parse_float(cols[16])
        d.ocean_linked    = cols[17] == "-" ? dh_label_t[] :
                            parse.(dh_label_t, split(cols[17], ','))
        deps[i] = d
    end

    return OracleDH(W, H, label, flowdirs, deps)
end

# Reproduces the C++ `main.cpp` ocean-labeling preamble (lines 57-73):
# label = NO_DEP, BucketFillFromEdges with ocean_level, then NaN -> OCEAN.
#
# Returns `(label, flowdirs, topo_clean)`. `topo_clean` is a copy of
# `topo` with NaN cells replaced by `-Inf` so the dephier priority queue
# orders them correctly (the C++ Array2D keeps the raw NoData numeric
# value like `-9999`; Phase 1's reader replaces it with NaN, which
# breaks heap ordering — substituting `-Inf` is equivalent and keeps the
# OCEAN wavefront popping first).
function prepare_label_and_flowdirs(topo_in::AbstractMatrix{<:AbstractFloat}, ocean_level)
    W, H = size(topo_in)
    topo = copy(topo_in)
    for i in eachindex(topo)
        if isnan(topo[i])
            topo[i] = -Inf
        end
    end
    label    = fill(NO_DEP,  W, H)
    flowdirs = fill(NO_FLOW, W, H)
    bucket_fill_from_edges!(topo, label, ocean_level, OCEAN)
    for i in eachindex(label)
        if isinf(topo[i]) && topo[i] < 0 || label[i] == OCEAN
            label[i] = OCEAN
        end
    end
    return label, flowdirs, topo
end

# Element-by-element comparison that treats NaN == NaN. Returns
# (matched::Bool, mismatch_count::Int, first_mismatch::Union{Nothing,Tuple}).
function diff_grid(name::AbstractString, jl, c)
    if size(jl) != size(c)
        return (false, length(jl), (; field = name, kind = :size, jl = size(jl), c = size(c)))
    end
    mismatches = 0
    first = nothing
    for i in eachindex(jl)
        if !isequal(jl[i], c[i])
            mismatches += 1
            if first === nothing
                first = (; field = name, idx = i, jl = jl[i], c = c[i])
            end
        end
    end
    return (mismatches == 0, mismatches, first)
end

# Compare two Depressions field by field. Returns Vector of
# (field, jl_value, c_value) for any mismatches.
function diff_depression(jl::Depression, c::Depression)
    mismatches = Tuple{Symbol, Any, Any}[]
    for f in (:pit_cell, :out_cell, :parent, :odep, :geolink,
              :pit_elev, :out_elev, :lchild, :rchild, :ocean_parent,
              :ocean_linked, :dep_label, :cell_count,
              :dep_vol, :water_vol, :total_elevation)
        jv = getfield(jl, f)
        cv = getfield(c,  f)
        if !isequal(jv, cv)
            push!(mismatches, (f, jv, cv))
        end
    end
    return mismatches
end
