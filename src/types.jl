# Mirrors the type definitions from
# submodules/dephier/include/dephier/dephier.hpp.
#
# Field names and types are kept verbatim to make the line-by-line port
# in dephier.jl read like the C++. `flat_c_idx` indices are stored 1-based
# in Julia (sentinel `NO_VALUE` = typemax(UInt32) is well above any
# realistic grid size, so 0-based↔1-based has no collision).

const dh_label_t = UInt32
const flat_c_idx = UInt32

mutable struct Depression{T}
    pit_cell::flat_c_idx
    out_cell::flat_c_idx
    parent::dh_label_t
    odep::dh_label_t
    geolink::dh_label_t
    pit_elev::T
    out_elev::T
    lchild::dh_label_t
    rchild::dh_label_t
    ocean_parent::Bool
    ocean_linked::Vector{dh_label_t}
    dep_label::dh_label_t
    cell_count::UInt32
    dep_vol::Float64
    water_vol::Float64
    total_elevation::Float64
end

# Default constructor matches the C++ in-class initializers.
function Depression{T}() where {T}
    return Depression{T}(
        NO_VALUE,
        NO_VALUE,
        NO_PARENT,
        NO_VALUE,
        NO_VALUE,
        T(Inf),
        T(Inf),
        NO_VALUE,
        NO_VALUE,
        false,
        dh_label_t[],
        dh_label_t(0),
        UInt32(0),
        0.0,
        0.0,
        0.0,
    )
end

const DepressionHierarchy{T} = Vector{Depression{T}}

# Outlet stores a candidate connection between two depressions. C++'s
# constructor swaps depa/depb so depa<=depb (preferred ordering for
# hashing); we do the same.
#
# We define the inner constructor explicitly with concrete field types
# so the convenience outer below doesn't recurse — Julia's auto-inner is
# `Outlet{T}(::Any, ::Any, ::Any, ::Any)`, which is *less* specific than
# `Outlet{T}(::Integer, ::Integer, ::Integer, ::T)`, so without an
# explicit typed inner the outer dispatches to itself forever.
mutable struct Outlet{T}
    depa::dh_label_t
    depb::dh_label_t
    out_cell::flat_c_idx
    out_elev::T

    function Outlet{T}(depa::dh_label_t, depb::dh_label_t, out_cell::flat_c_idx, out_elev::T) where {T}
        return new{T}(depa, depb, out_cell, out_elev)
    end
end

function Outlet{T}(depa::Integer, depb::Integer, out_cell::Integer, out_elev::T) where {T}
    if depa > depb
        depa, depb = depb, depa
    end
    return Outlet{T}(dh_label_t(depa), dh_label_t(depb), flat_c_idx(out_cell), out_elev)
end

# In C++ this is `OutletLink` — a separate struct used as a key in
# `unordered_map<OutletLink, Outlet, OutletHash>`. Julia's `Dict` keyed
# on a tuple is identical in semantics, so we drop the wrapper.
const OutletLink = Tuple{dh_label_t, dh_label_t}
