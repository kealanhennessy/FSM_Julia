# Phase 4 randomized property tests — direct ports of the random
# TEST_CASEs in `unittests/fsm_tests.cpp`. Where the upstream uses
# perlin-on-the-fly with std::mt19937_64, we use a fixed batch of
# pre-generated TIFFs under test/test_cases/random/ (one terrain per
# iteration). That keeps the tests deterministic without needing a
# Julia port of perlin.
#
# Properties exercised:
#   * `MoveWaterIntoPits Repeatedly` — seeding wtd with the previous
#     run's pit-cell volumes reproduces the per-depression accumulation.
#   * `Randomized Heavy Flooding vs Priority-Flood` — FSM with very
#     large initial wtd produces the same flooded surface as Zhou2016
#     Priority-Flood applied to the topography alone.
#   * `Randomized Testing of Repeated FSM` — applying FSM twice in a
#     row leaves wtd unchanged after the second run.
#   * `Randomized Testing of Incremental FSM vs Big Dump` — one bulk
#     fill produces the same wtd as ten incremental 0.1 fills.
#   * `RandomizedMassConservation` — total water (sum(wtd) + ocean
#     depression's water_vol) is conserved against the input.

import FillSpillMerge:
    OCEAN, NO_DEP, NO_VALUE, NO_FLOW,
    Depression, dh_label_t,
    get_depression_hierarchy,
    fill_spill_merge!,
    move_water_into_pits!,
    priority_flood_zhou2016!

# All terrain files live here. Filename prefix encodes (size class) ×
# (elev kind) — see tools/gen_random_terrains.cpp.
const RANDOM_TERRAINS_DIR = joinpath(@__DIR__, "test_cases", "random")

# Lazy lazily-loaded inventory of terrain paths. We iterate over the
# smaller / integer set by default for speed; the heavy-flooding test
# uses the float set since the C++ uses non-truncated perlin there.
function _list_terrains(pattern::AbstractString)
    isdir(RANDOM_TERRAINS_DIR) || return String[]
    return sort([
        joinpath(RANDOM_TERRAINS_DIR, f)
        for f in readdir(RANDOM_TERRAINS_DIR)
        if startswith(f, pattern)
    ])
end

# Load a TIFF as a (W, H) Float64 matrix, preserving the -1.0 ocean
# ring. We deliberately do NOT call read_tif (which replaces NoData with
# NaN) — the generator's NoData is set to -9999, but no cell in the
# data actually carries that value; the ocean ring is numerically -1.0
# and we want to keep it that way for the algorithm.
function _load_random_terrain(path::AbstractString)
    return ArchGDAL.read(path) do dataset
        band = ArchGDAL.getband(dataset, 1)
        arr  = ArchGDAL.read(band)
        Matrix{Float64}(arr)
    end
end

# Build (label, flowdirs) for a random-batch terrain: outer ring is
# OCEAN, everything else is NO_DEP / NO_FLOW. Mirrors C++
# `label.setEdges(OCEAN)`.
function _init_label_flowdirs(W::Integer, H::Integer)
    label    = fill(NO_DEP,  W, H)
    flowdirs = fill(NO_FLOW, W, H)
    set_edges_ocean!(label)
    return label, flowdirs
end


@testset "C++ port: MoveWaterIntoPits Repeatedly" begin
    # Mirrors `TEST_CASE("MoveWaterIntoPits Repeatedly")`. Property:
    # accumulating wtd=1 over the whole grid, then re-seeding wtd with
    # the resulting per-depression volumes at pit cells, then running
    # MoveWaterIntoPits a second time, gives the same per-depression
    # water_vol.
    paths = vcat(_list_terrains("small_int_"), _list_terrains("large_int_"))
    @test !isempty(paths)
    for path in paths
        topo = _load_random_terrain(path)
        W, H = size(topo)

        # First pass: get the hierarchy, accumulate wtd=1 into pits.
        label, flowdirs = _init_label_flowdirs(W, H)
        deps1 = get_depression_hierarchy(topo, label, flowdirs)
        wtd = fill(1.0, W, H)
        move_water_into_pits!(topo, label, flowdirs, deps1, wtd)

        # Second pass: fresh hierarchy on a freshly-initialised label,
        # then seed wtd at pit cells with the per-dep water_vol from
        # the first pass.
        label, flowdirs = _init_label_flowdirs(W, H)
        deps2 = get_depression_hierarchy(topo, label, flowdirs)
        wtd = zeros(W, H)
        for dep in deps1
            if dep.water_vol > 0 && dep.pit_cell != NO_VALUE
                wtd[dep.pit_cell] = dep.water_vol
            end
        end
        move_water_into_pits!(topo, label, flowdirs, deps2, wtd)

        # Per-depression water_vol should match (skipping OCEAN at
        # index 1).
        for i in 2:length(deps1)
            @test deps1[i].water_vol == deps2[i].water_vol
        end
    end
end


@testset "C++ port: Randomized Testing of Repeated FSM" begin
    # Mirrors `TEST_CASE("Randomized Testing of Repeated FSM")`.
    # Property: FSM is idempotent on its own output — running it a
    # second time leaves wtd unchanged.
    paths = vcat(_list_terrains("small_int_"), _list_terrains("large_int_"))
    @test !isempty(paths)
    for path in paths
        topo_orig = _load_random_terrain(path)
        W, H = size(topo_orig)

        wtd = fill(1.0, W, H)

        function do_fsm!()
            topo = copy(topo_orig)  # FSM doesn't mutate topo, but be safe
            label, flowdirs = _init_label_flowdirs(W, H)
            deps = get_depression_hierarchy(topo, label, flowdirs)
            fill_spill_merge!(topo, label, flowdirs, deps, wtd)
        end

        do_fsm!()
        first_wtd = copy(wtd)
        do_fsm!()

        @test maximum(abs.(first_wtd .- wtd)) < 1e-6
    end
end


@testset "C++ port: Randomized Testing of Incremental FSM vs Big Dump" begin
    # Mirrors `TEST_CASE("Randomized Testing of Incremental FSM vs Big Dump")`.
    # Property: one bulk fill (swl=1) gives the same wtd as ten
    # incremental fills of 0.1.
    paths = vcat(_list_terrains("small_int_"), _list_terrains("large_int_"))
    @test !isempty(paths)
    for path in paths
        topo_orig = _load_random_terrain(path)
        W, H = size(topo_orig)

        # Big dump.
        topo = copy(topo_orig)
        label, flowdirs = _init_label_flowdirs(W, H)
        deps = get_depression_hierarchy(topo, label, flowdirs)
        wtd_big = fill(1.0, W, H)
        fill_spill_merge!(topo, label, flowdirs, deps, wtd_big)

        # Ten incremental fills of 0.1.
        wtd_inc = zeros(W, H)
        for _step in 1:10
            wtd_inc .+= 0.1
            topo = copy(topo_orig)
            label, flowdirs = _init_label_flowdirs(W, H)
            deps = get_depression_hierarchy(topo, label, flowdirs)
            fill_spill_merge!(topo, label, flowdirs, deps, wtd_inc)
        end

        @test maximum(abs.(wtd_big .- wtd_inc)) < 1e-6
    end
end


@testset "C++ port: RandomizedMassConservation" begin
    # Mirrors `TEST_CASE("RandomizedMassConservation")`. Property: water
    # is conserved. Total above-ground water in the interior + water
    # routed to the OCEAN depression equals the input swl × cells. We
    # use a fixed sequence of swl values derived from the path's seed
    # for reproducibility.
    paths = vcat(_list_terrains("small_int_"), _list_terrains("large_int_"))
    @test !isempty(paths)
    for (i, path) in enumerate(paths)
        topo = _load_random_terrain(path)
        W, H = size(topo)

        # Deterministic surface_water_amount in (0, 1] keyed by index.
        # Matches the C++'s uniform_real_distribution(0,1) draws in
        # spirit (not in value) — we only need a varied set.
        swl = (i + 0.5) / (length(paths) + 1)

        label, flowdirs = _init_label_flowdirs(W, H)
        deps = get_depression_hierarchy(topo, label, flowdirs)
        wtd = fill(swl, W, H)

        fill_spill_merge!(topo, label, flowdirs, deps, wtd)

        total = sum(wtd) + deps[1].water_vol  # deps[1] is OCEAN
        expected = swl * W * H
        # doctest::Approx default is 1e-5 relative; we use a loose
        # absolute tolerance scaled by the input total.
        @test abs(total - expected) < 1e-6 * max(1.0, abs(expected))
    end
end


@testset "C++ port: Randomized Heavy Flooding vs Priority-Flood" begin
    # Mirrors `TEST_CASE("Randomized Heavy Flooding vs Priority-Flood")`.
    # Property: with abundant initial water (wtd = 100), every
    # depression is filled to its brim. The resulting hydrologic surface
    # (topo + wtd) should match what Zhou2016 Priority-Flood produces
    # when applied to the topography alone (since Priority-Flood floods
    # everything as if water were abundant).
    #
    # Uses float terrains (not integer-truncated) per the C++.
    paths = vcat(_list_terrains("small_float_"), _list_terrains("large_float_"))
    @test !isempty(paths)
    for path in paths
        topo_orig = _load_random_terrain(path)
        W, H = size(topo_orig)

        # FSM heavy flood.
        topo = copy(topo_orig)
        label, flowdirs = _init_label_flowdirs(W, H)
        deps = get_depression_hierarchy(topo, label, flowdirs)
        wtd = fill(100.0, W, H)
        fill_spill_merge!(topo, label, flowdirs, deps, wtd)
        fsm_surface = topo .+ wtd

        # Priority-Flood reference. Operates in-place; copy first.
        pf_surface = copy(topo_orig)
        priority_flood_zhou2016!(pf_surface)

        # The C++ tolerance is 1e-6; PF and FSM both leave ocean cells
        # at their original (numerical) elevation, so the comparison
        # spans the whole grid.
        @test maximum(abs.(fsm_surface .- pf_surface)) < 1e-6
    end
end
