# Phase 4 unit tests — direct ports of unittests/fsm_tests.cpp from the
# upstream Barnes2020-FillSpillMerge repo. Each Julia testset corresponds
# 1-to-1 with a C++ TEST_CASE; SUBCASEs map to nested @testsets.
#
# Conventions for the test fixtures:
# - C++ matrix literals are written row-by-row of the visual grid. In
#   Julia we'd naturally write the same syntax, but `M[i, j]` in Julia
#   gives "row i, col j" while our project convention is `M[x, y]` =
#   "col x, row y" (mirroring C++ Array2D). So we write each fixture as
#   a Julia matrix literal in visual order, then call `visual_grid` to
#   transpose into the (W, H) layout the algorithm expects. The
#   resulting `topo[x, y]` matches the C++ source's `topo(x, y)`.
# - The C++ uses 0-based `xyToI(x, y)`; Julia is 1-based, so coordinates
#   from the C++ source are bumped by 1 before being passed to
#   `LinearIndices`. A helper `xy_to_flat` does the bump.
# - `dh_label_t` is `UInt32`; `flowdir_t` is `Int8`; both match the
#   richdem source.

import FillSpillMerge:
    fp_eq, fp_le, fp_ge,
    D8X, D8Y, D8_INVERSE,
    OCEAN, NO_DEP, NO_VALUE, NO_PARENT, NO_FLOW,
    Depression, dh_label_t, flat_c_idx,
    DepressionHierarchy,
    get_depression_hierarchy,
    fill_spill_merge!,
    depression_volume,
    move_water_into_pits!,
    backfill_depression!,
    determine_water_level!,
    fill_depressions!

# Transpose a visually-laid-out matrix literal into the algorithm's
# (W, H) `M[x, y]` layout.
visual_grid(m::AbstractMatrix) = permutedims(m)

# Convert C++ 0-based (x, y) to Julia 1-based flat linear index for an
# array of the given (W, H) size.
xy_to_flat(W::Integer, H::Integer, x::Integer, y::Integer) =
    Int(LinearIndices((W, H))[x + 1, y + 1])

# Set the outer ring of a label grid to OCEAN. Mirrors C++
# `label.setEdges(OCEAN)`.
function set_edges_ocean!(label::AbstractMatrix)
    W, H = size(label)
    @inbounds for x in 1:W
        label[x, 1] = OCEAN
        label[x, H] = OCEAN
    end
    @inbounds for y in 1:H
        label[1, y] = OCEAN
        label[W, y] = OCEAN
    end
    return label
end


@testset "C++ port: Depression volume" begin
    # Mirrors `TEST_CASE("Depression volume")`.
    @test depression_volume(2, 5, 10) == 0
    @test depression_volume(3, 5, 10) == 5
    @test depression_volume(4, 5, 10) == 10
end


@testset "C++ port: Determine water level" begin
    # Mirrors `TEST_CASE("Determine water level")`. The C++ takes
    # `sill_wtd` by reference; we pass `(wtd, cx, cy)` and read out the
    # mutated cell afterwards.
    @testset "depression volume exactly equals water volume" begin
        wtd = fill(0.0, 1, 1)
        wtd[1, 1] = -2.0
        water_level = determine_water_level!(wtd, 1, 1, 10.0, 4.0, 5, 10.0)
        @test wtd[1, 1] == -2.0
        @test water_level == 4.0
    end
    @testset "water volume less than depression volume" begin
        wtd = fill(0.0, 1, 1)
        wtd[1, 1] = -2.0
        water_level = determine_water_level!(wtd, 1, 1, 8.0, 4.0, 5, 10.0)
        @test wtd[1, 1] == -2.0
        @test water_level == 18 / 5
    end
    @testset "water volume greater than depression volume" begin
        wtd = fill(0.0, 1, 1)
        wtd[1, 1] = -2.0
        water_level = determine_water_level!(wtd, 1, 1, 12.0, 4.0, 5, 10.0)
        @test wtd[1, 1] == 0.0
        @test water_level == 4.0
    end
end


@testset "C++ port: MoveWaterIntoPits 1" begin
    # Mirrors `TEST_CASE("MoveWaterIntoPits 1")`. 10x10 grid; pit at
    # (7, 7) (C++ 0-based: 7, 7).
    topo = visual_grid([
        -9 -9 -9 -9 -9 -9 -9 -9 -9 -9;
        -9  9  9  9  9  9  9  9  9 -9;
        -9  9  8  8  8  8  7  6  9 -9;
        -9  9  8  7  8  7  6  5  9 -9;
        -9  9  8  7  8  6  5  4  9 -9;
        -9  9  8  8  8  5  4  3  9 -9;
        -9  9  7  6  5  4  3  2  9 -9;
        -9  9  7  6  5  4  3  1  9 -9;
        -9  9  9  9  9  9  9  9  9 -9;
        -9 -9 -9 -9 -9 -9 -9 -9 -9 -9;
    ]) .|> Float64

    W, H = size(topo)
    label    = fill(NO_DEP,  W, H)
    set_edges_ocean!(label)
    flowdirs = fill(NO_FLOW, W, H)
    wtd      = fill(0.0,     W, H)

    deps = get_depression_hierarchy(topo, label, flowdirs)

    label_good = visual_grid([
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
        0 0 2 2 1 1 1 1 0 0;
        0 0 2 2 1 1 1 1 0 0;
        0 0 2 2 1 1 1 1 0 0;
        0 0 1 1 1 1 1 1 0 0;
        0 0 1 1 1 1 1 1 0 0;
        0 0 1 1 1 1 1 1 0 0;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
    ]) .|> dh_label_t

    flowdirs_good = visual_grid([
        0 0 0 0 0 0 0 0 0 0;
        0 8 4 4 4 4 4 4 6 0;
        0 8 6 7 6 6 6 7 6 0;
        0 8 6 7 6 6 6 7 6 0;
        0 8 5 0 6 6 6 7 6 0;
        0 8 6 6 6 6 6 7 6 0;
        0 8 6 5 6 5 6 7 6 0;
        0 8 5 4 5 4 5 0 6 0;
        0 6 6 6 6 6 6 6 6 0;
        0 0 0 0 0 0 0 0 0 0;
    ]) .|> Int8

    @test label    == label_good
    @test flowdirs == flowdirs_good

    fill!(wtd, 1.0)
    move_water_into_pits!(topo, label, flowdirs, deps, wtd)

    @test all(==(0.0), wtd)

    @test deps[1].water_vol == 64.0    # OCEAN (label 0)
    @test deps[2].water_vol == 30.0    # label 1
    @test deps[3].water_vol ==  6.0    # label 2

    @test deps[1].parent == NO_VALUE   # OCEAN
    @test deps[2].parent == 3
    @test deps[3].parent == 3
    @test deps[4].parent == 0

    @test isnan(deps[1].dep_vol)       # OCEAN: 0 * Inf = NaN, same as C++
    @test deps[2].dep_vol ==  73.0
    @test deps[3].dep_vol ==   2.0
    @test deps[4].dep_vol == 111.0
end


@testset "C++ port: Backfill Depression" begin
    # Mirrors `TEST_CASE("Backfill Depression")`. Fixed cells_affected
    # list, expected wtd_good.
    topo = visual_grid([
        -9 -9 -9 -9 -9 -9 -9 -9 -9 -9;
        -9  6  6  6  6  6  6  6  6 -9;
        -9  6  1  6  1  1  6  1  6 -9;
        -9  6  1  6  1  3  6  1  6 -9;
        -9  6  1  6  2  1  4  1  6 -9;
        -9  6  1  6  1  1  6  1  6 -9;
        -9  6  1  6  6  6  6  1  6 -9;
        -9  6  1  1  1  1  1  1  6 -9;
        -9  6  6  6  6  6  6  6  6 -9;
        -9 -9 -9 -9 -9 -9 -9 -9 -9 -9;
    ]) .|> Float64
    W, H = size(topo)

    wtd = zeros(W, H)

    cells_affected = [
        xy_to_flat(W, H, 4, 2), xy_to_flat(W, H, 5, 2),
        xy_to_flat(W, H, 4, 3), xy_to_flat(W, H, 5, 3),
        xy_to_flat(W, H, 4, 4), xy_to_flat(W, H, 5, 4),
        xy_to_flat(W, H, 4, 5), xy_to_flat(W, H, 5, 5),
    ]

    backfill_depression!(4.0, topo, wtd, cells_affected)

    wtd_good = visual_grid([
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 3 3 0 0 0 0;
        0 0 0 0 3 1 0 0 0 0;
        0 0 0 0 2 3 0 0 0 0;
        0 0 0 0 3 3 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
    ]) .|> Float64
    @test wtd == wtd_good
end


@testset "C++ port: FillDepressions" begin
    # Mirrors `TEST_CASE("FillDepressions")`. 10x10 grid; 5 SUBCASEs.
    topo = visual_grid([
        -9 -9 -9 -9 -9 -9 -9 -9 -9 -9;
        -9  6  6  6  6  6  6  6  6 -9;
        -9  6  1  6  1  1  6  1  6 -9;
        -9  6  1  6  3  3  6  1  6 -9;
        -9  6  1  6  2  1  4  1  6 -9;
        -9  6  1  6  1  1  6  1  6 -9;
        -9  6  1  6  6  6  6  1  6 -9;
        -9  6  1  1  1  1  1  1  6 -9;
        -9  6  6  6  6  6  6  6  6 -9;
        -9 -9 -9 -9 -9 -9 -9 -9 -9 -9;
    ]) .|> Float64
    label = visual_grid([
        0 0 0 0 0 0 0 0 0 0;
        0 1 1 2 2 2 2 1 1 0;
        0 1 1 2 2 2 2 1 1 0;
        0 1 1 2 2 2 2 1 1 0;
        0 1 1 2 2 2 2 1 1 0;
        0 1 1 2 2 2 2 1 1 0;
        0 1 1 2 2 2 2 1 1 0;
        0 1 1 1 1 1 1 1 1 0;
        0 1 1 1 1 1 1 1 1 0;
        0 0 0 0 0 0 0 0 0 0;
    ]) .|> dh_label_t
    W, H = size(topo)

    pit_cell = xy_to_flat(W, H, 4, 2)
    out_cell = xy_to_flat(W, H, 4, 3)

    @test topo[pit_cell] == 1.0

    dep_labels = Set{dh_label_t}([dh_label_t(2)])

    @testset "no water to add" begin
        wtd = zeros(W, H)
        wtd[5, 4] = -0.5  # C++ wtd(4, 3) -> Julia 1-based [5, 4]
        wtd_good = copy(wtd)
        fill_depressions!(pit_cell, out_cell, dep_labels, 0.0, topo, label, wtd)
        @test wtd == wtd_good
    end

    @testset "standard case" begin
        wtd = zeros(W, H)
        fill_depressions!(pit_cell, out_cell, dep_labels, 3.0, topo, label, wtd)
        W_ = 1.5
        wtd_good = visual_grid([
            0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 W_ W_ 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
        ]) .|> Float64
        @test maximum(abs.(wtd .- wtd_good)) < 1e-6
    end

    @testset "sill cell absorbs some water" begin
        wtd = zeros(W, H)
        wtd[5, 4] = -1.0   # C++ wtd(4, 3)
        fill_depressions!(pit_cell, out_cell, dep_labels, 5.0, topo, label, wtd)
        W_ = 2.0
        wtd_good = visual_grid([
            0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 W_ W_ 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
            0 0 0 0 0  0 0 0 0 0;
        ]) .|> Float64
        @test maximum(abs.(wtd .- wtd_good)) < 1e-6
    end

    @testset "passes over a saddle" begin
        wtd = zeros(W, H)
        # C++ uses out_cell = topo.xyToI(6, 4) for this subcase.
        out_cell_2 = xy_to_flat(W, H, 6, 4)
        fill_depressions!(pit_cell, out_cell_2, dep_labels, 19.0, topo, label, wtd)
        wtd_good = visual_grid([
            0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 3 3 0 0 0 0;
            0 0 0 0 1 1 0 0 0 0;
            0 0 0 0 2 3 0 0 0 0;
            0 0 0 0 3 3 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0;
        ]) .|> Float64
        @test maximum(abs.(wtd .- wtd_good)) < 1e-6
    end

    @testset "passes over a saddle and sill absorbs some" begin
        wtd = zeros(W, H)
        out_cell_2 = xy_to_flat(W, H, 6, 4)
        wtd[7, 5] = -1.0   # C++ wtd(6, 4) -> Julia [7, 5]
        fill_depressions!(pit_cell, out_cell_2, dep_labels, 19.5, topo, label, wtd)
        wtd_good = visual_grid([
            0 0 0 0 0 0    0 0 0 0;
            0 0 0 0 0 0    0 0 0 0;
            0 0 0 0 3 3    0 0 0 0;
            0 0 0 0 1 1    0 0 0 0;
            0 0 0 0 2 3 -0.5 0 0 0;
            0 0 0 0 3 3    0 0 0 0;
            0 0 0 0 0 0    0 0 0 0;
            0 0 0 0 0 0    0 0 0 0;
            0 0 0 0 0 0    0 0 0 0;
            0 0 0 0 0 0    0 0 0 0;
        ]) .|> Float64
        @test maximum(abs.(wtd .- wtd_good)) < 1e-6
    end
end


@testset "C++ port: PQ Issue" begin
    # Mirrors `TEST_CASE("PQ Issue")`. Regression test: just verifies
    # FSM doesn't crash / throw on this specific layout.
    topo = visual_grid([
        -9 -9 -9 -9 -9 -9 -9  9 -9;
        -9  5  5  5  5  5  5  5 -9;
        -9  5  1  1  5  1  1  5 -9;
        -9  5  1  1  5  1  1  5 -9;
        -9  2  1  1  3  1  1  5 -9;
        -9  5  1  1  5  1  1  5 -9;
        -9  5  5  5  5  5  5  5 -9;
        -9 -9 -9 -9 -9 -9 -9  9 -9;
    ]) .|> Float64
    W, H = size(topo)

    label    = fill(NO_DEP,  W, H)
    set_edges_ocean!(label)
    flowdirs = fill(NO_FLOW, W, H)

    deps = get_depression_hierarchy(topo, label, flowdirs)
    wtd  = zeros(W, H)
    wtd[6, 6] = 17.0   # C++ wtd(5, 5) -> Julia [6, 6]

    # The original test has no output assertion; it's a regression
    # against a bug that used to throw. We assert non-throw.
    @test (fill_spill_merge!(topo, label, flowdirs, deps, wtd); true)
end


@testset "C++ port: PQ Issue 2" begin
    # Mirrors `TEST_CASE("PQ Issue 2")`.
    topo = visual_grid([
        -9 -9 -9 -9 -9;
        -9  0  0  0 -9;
        -9  0 -8  0 -9;
        -9  0 -8  0 -9;
        -9  0  0  0 -9;
        -9 -9 -9 -9 -9;
    ]) .|> Float64
    labels = visual_grid([
        0 0 0 0 0;
        0 0 0 0 0;
        0 0 1 0 0;
        0 0 1 0 0;
        0 0 0 0 0;
        0 0 0 0 0;
    ]) .|> dh_label_t
    W, H = size(topo)

    pit_cell = xy_to_flat(W, H, 2, 2)
    out_cell = xy_to_flat(W, H, 3, 3)

    dep_labels = Set{dh_label_t}([dh_label_t(1)])

    wtd = zeros(W, H)
    @test (fill_depressions!(pit_cell, out_cell, dep_labels, 2.0, topo, labels, wtd); true)
end
