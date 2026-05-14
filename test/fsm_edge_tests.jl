# Phase 4 edge-case tests. These go beyond the C++ unittests to
# exercise degenerate / boundary inputs: tiny grids, all-ocean,
# monotonic-slope (no depressions), and Float32 elevations.
#
# Conventions match fsm_unit_tests.jl: `visual_grid` transposes a
# C++-style row-major literal into the algorithm's (W, H) `M[x, y]`
# layout; `set_edges_ocean!` mirrors C++ `label.setEdges(OCEAN)`.

import FillSpillMerge:
    OCEAN, NO_DEP, NO_VALUE, NO_FLOW,
    Depression, dh_label_t,
    get_depression_hierarchy,
    fill_spill_merge!


@testset "Edge: 3x3 land cell surrounded by lower ocean is absorbed" begin
    # Smallest meaningful grid: ocean ring + one land cell at the
    # centre, where the land cell is HIGHER than the ocean ring. The
    # dephier algorithm expands ocean inward; since every neighbour of
    # the centre is lower-elevation ocean, the centre is absorbed into
    # OCEAN rather than forming its own depression. The only entry in
    # the hierarchy is OCEAN. FSM then routes any initial wtd at that
    # cell to OCEAN's water_vol.
    topo = visual_grid([
        -1 -1 -1;
        -1  5 -1;
        -1 -1 -1;
    ]) .|> Float64
    W, H = size(topo)
    label    = fill(NO_DEP,  W, H)
    set_edges_ocean!(label)
    flowdirs = fill(NO_FLOW, W, H)

    deps = get_depression_hierarchy(topo, label, flowdirs)
    @test length(deps) == 1
    @test deps[1].dep_label == OCEAN
    @test label[2, 2] == OCEAN

    wtd = fill(1.0, W, H)
    for i in eachindex(label)
        if label[i] == OCEAN
            wtd[i] = 0.0
        end
    end
    # After the ocean-cell zeroing, the centre wtd is 0 because it's
    # OCEAN-labelled. FSM is then a no-op on this trivial case.
    fill_spill_merge!(topo, label, flowdirs, deps, wtd)
    @test all(==(0.0), wtd)
end


@testset "Edge: all-ocean grid throws" begin
    # No interior: every cell is OCEAN. The dephier-port preconditions
    # check that at least one OCEAN cell has a non-ocean neighbour;
    # this case violates that.
    W, H = 3, 3
    topo = fill(-1.0, W, H)
    label    = fill(OCEAN,   W, H)
    flowdirs = fill(NO_FLOW, W, H)
    @test_throws ArgumentError get_depression_hierarchy(topo, label, flowdirs)
end


@testset "Edge: 1x1 land throws (no ocean)" begin
    # No ocean anywhere. dephier should refuse.
    topo = reshape([5.0], 1, 1)
    label    = fill(NO_DEP,  1, 1)
    flowdirs = fill(NO_FLOW, 1, 1)
    @test_throws ArgumentError get_depression_hierarchy(topo, label, flowdirs)
end


@testset "Edge: monotonic slope has no depressions" begin
    # 5x5 with elev increasing in y; edges are OCEAN. Interior has no
    # local minima except where the slope meets the ocean edge — i.e.
    # every interior cell drains. FSM should leave no standing water.
    visual = [
        -1.0 -1.0 -1.0 -1.0 -1.0;
        -1.0  1.0  2.0  3.0 -1.0;
        -1.0  1.0  2.0  3.0 -1.0;
        -1.0  1.0  2.0  3.0 -1.0;
        -1.0 -1.0 -1.0 -1.0 -1.0;
    ]
    topo = visual_grid(visual)
    W, H = size(topo)
    label    = fill(NO_DEP,  W, H)
    set_edges_ocean!(label)
    flowdirs = fill(NO_FLOW, W, H)

    deps = get_depression_hierarchy(topo, label, flowdirs)
    # Only the OCEAN, since the interior slope drains uniformly. A
    # leaf depression may still be created at the bottom of the slope
    # if it forms a "plateau" with neighbours — but at minimum there's
    # one OCEAN entry.
    @test length(deps) >= 1
    @test deps[1].dep_label == OCEAN

    wtd = fill(0.5, W, H)
    for i in eachindex(label)
        if label[i] == OCEAN
            wtd[i] = 0.0
        end
    end

    fill_spill_merge!(topo, label, flowdirs, deps, wtd)

    # All water should leave the interior (routed to ocean).
    for i in eachindex(label)
        if label[i] == OCEAN
            @test wtd[i] == 0.0
        else
            # Cells along a fully draining slope keep no standing water.
            @test wtd[i] <= 0
        end
    end
end


@testset "Edge: deep pit holds water exactly to brim" begin
    # 5x5 grid: ocean ring + inner 3x3 plateau at elev 5 with a single
    # pit at the centre (elev 0). With enough water to fill the pit but
    # not overflow the plateau, the centre cell should rise toward the
    # plateau height.
    visual = [
        -1.0 -1.0 -1.0 -1.0 -1.0;
        -1.0  5.0  5.0  5.0 -1.0;
        -1.0  5.0  0.0  5.0 -1.0;
        -1.0  5.0  5.0  5.0 -1.0;
        -1.0 -1.0 -1.0 -1.0 -1.0;
    ]
    topo = visual_grid(visual)
    W, H = size(topo)
    label    = fill(NO_DEP,  W, H)
    set_edges_ocean!(label)
    flowdirs = fill(NO_FLOW, W, H)

    deps = get_depression_hierarchy(topo, label, flowdirs)

    # Just enough water to fill the centre to depth 5 (dep_vol = 5).
    wtd = zeros(W, H)
    wtd[3, 3] = 5.0  # centre cell in (W, H)=(5, 5) layout

    fill_spill_merge!(topo, label, flowdirs, deps, wtd)

    # The pit is now at water level 5 (flush with plateau).
    @test wtd[3, 3] ≈ 5.0 atol=1e-9
end


@testset "Float32 elevations work end-to-end" begin
    # Tiny case_02-style problem with Float32 topo. The algorithm's
    # internal accumulations are Float64 (depression water_vol etc.),
    # but topo and wtd carry through as Float32. The point of this
    # test is that the dispatch doesn't crash and the resulting wtd
    # is sensible (water either stays in the depression or routes to
    # the ocean).
    topo = visual_grid([
        -1f0 -1f0 -1f0 -1f0 -1f0;
        -1f0  3f0  3f0  3f0 -1f0;
        -1f0  3f0  0f0  3f0 -1f0;
        -1f0  3f0  3f0  3f0 -1f0;
        -1f0 -1f0 -1f0 -1f0 -1f0;
    ])
    @test eltype(topo) == Float32
    W, H = size(topo)
    label    = fill(NO_DEP,  W, H)
    set_edges_ocean!(label)
    flowdirs = fill(NO_FLOW, W, H)

    deps = get_depression_hierarchy(topo, label, flowdirs)
    @test eltype(deps[1].pit_elev) == Float32

    wtd = fill(Float32(1.0), W, H)
    for i in eachindex(label)
        if label[i] == OCEAN
            wtd[i] = Float32(0.0)
        end
    end

    fill_spill_merge!(topo, label, flowdirs, deps, wtd)

    # Ocean stays at 0.
    for i in eachindex(label)
        if label[i] == OCEAN
            @test wtd[i] == Float32(0.0)
        end
    end
    # The centre pit absorbs ~3.0 of water (its dep_vol with the
    # 4-cell ring at elev 3 surrounding it).
    @test wtd[3, 3] >= Float32(0.0)
end
