# Mirrors submodules/dephier/include/dephier/DisjointDenseIntSet.hpp.
#
# Two main entry points used by dephier:
#   - `find_set!`  (path compression by mutation; mirrors `findSet`)
#   - `merge_a_into_b!` (preserves parenthood; A becomes child of B)
#
# Indices are 0-based to match the C++ — depression labels (`dh_label_t`)
# already start at 0 (the OCEAN), so keeping the underlying array
# 0-indexed via 1-based Julia storage `parent[i+1]` would just add noise.
# Instead we offset internally: store `parent[1+i]` and expose i.

mutable struct DisjointDenseIntSet
    rank::Vector{UInt32}
    parent::Vector{UInt32}
end

function DisjointDenseIntSet()
    return DisjointDenseIntSet(UInt32[], UInt32[])
end

# Pre-size with N initial sets, each its own parent.
function DisjointDenseIntSet(N::Integer)
    rank   = zeros(UInt32, N)
    parent = UInt32.(0:N-1)
    return DisjointDenseIntSet(rank, parent)
end

# Grow the data structure so that index `n` is valid; new sets are
# their own parent with rank 0.
function check_size!(s::DisjointDenseIntSet, n::Integer)
    needed = Int(n) + 1
    if needed <= length(s.rank)
        return
    end
    old_size = length(s.rank)
    resize!(s.rank,   needed); s.rank[old_size+1:end]   .= 0
    resize!(s.parent, needed)
    for i in old_size:needed-1
        s.parent[i+1] = i
    end
    return
end

# Explicitly create a set (and any intermediates).
make_set!(s::DisjointDenseIntSet, n::Integer) = check_size!(s, n)

max_element(s::DisjointDenseIntSet) = length(s.rank) - 1

# Returns the representative of the set containing n, applying path
# compression so subsequent queries are O(1).
function find_set!(s::DisjointDenseIntSet, n::Integer)
    if n + 1 > length(s.parent)
        throw(ArgumentError("find_set!($n) is outside the valid range [0,$(length(s.parent)))"))
    end
    p = s.parent[n+1]
    if p == n
        return UInt32(n)
    else
        root = find_set!(s, Int(p))
        s.parent[n+1] = root
        return root
    end
end

# Union by rank. Resulting parent is unpredictable.
function union_set!(s::DisjointDenseIntSet, a::Integer, b::Integer)
    roota = find_set!(s, a)
    rootb = find_set!(s, b)
    if roota == rootb
        return
    end
    if s.rank[roota+1] < s.rank[rootb+1]
        s.parent[roota+1] = rootb
    elseif s.rank[roota+1] > s.rank[rootb+1]
        s.parent[rootb+1] = roota
    else
        s.parent[rootb+1] = roota
        s.rank[roota+1] += 1
    end
    return
end

# Force A to become a child of B. Sacrifices balance for predictable
# parenthood — used when the dephier hierarchy needs the new metanode
# to be the parent of the two children we just merged.
function merge_a_into_b!(s::DisjointDenseIntSet, a::Integer, b::Integer)
    check_size!(s, a)
    check_size!(s, b)
    s.parent[a+1] = b
    if s.rank[a+1] == s.rank[b+1]
        s.rank[b+1] += 1
    elseif s.rank[a+1] > s.rank[b+1]
        s.rank[b+1] = s.rank[a+1] + 1
    end
    return
end

same_set(s::DisjointDenseIntSet, a::Integer, b::Integer) = find_set!(s, a) == find_set!(s, b)
