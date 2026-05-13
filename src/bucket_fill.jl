# Mirrors richdem/misc/misc_methods.hpp BucketFill / BucketFillFromEdges.
#
# Note: the check is exact equality (`!=` against `check_value`), not `<=`.
# For FSM's ocean detection, this means edge-bucket-fill alone labels only
# cells whose elevation is exactly `ocean_level`; the NoData -> OCEAN step
# in main.cpp covers everything else.

function bucket_fill!(check, set, check_value, set_value, seeds)
    W, H = size(check)
    while !isempty(seeds)
        x, y = pop!(seeds)
        if check[x, y] != check_value || set[x, y] == set_value
            continue
        end
        set[x, y] = set_value
        for n in 1:8
            nx, ny = x + D8X[n], y + D8Y[n]
            if !(1 <= nx <= W && 1 <= ny <= H)
                continue
            end
            if check[nx, ny] == check_value && set[nx, ny] != set_value
                push!(seeds, (nx, ny))
            end
        end
    end
    return set
end

function bucket_fill_from_edges!(check, set, check_value, set_value)
    W, H = size(check)
    size(set) == (W, H) || throw(DimensionMismatch("check and set must have the same dimensions"))
    seeds = Tuple{Int,Int}[]
    sizehint!(seeds, 2W + 2H)
    for y in 1:H
        push!(seeds, (1, y))
        push!(seeds, (W, y))
    end
    for x in 1:W
        push!(seeds, (x, 1))
        push!(seeds, (x, H))
    end
    return bucket_fill!(check, set, check_value, set_value, seeds)
end
