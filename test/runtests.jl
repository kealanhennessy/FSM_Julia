using Test
using FillSpillMerge
using ArchGDAL

import FillSpillMerge:
    fp_eq, fp_le, fp_ge,
    bucket_fill!, bucket_fill_from_edges!,
    D8X, D8Y, D8_INVERSE,
    OCEAN, NO_DEP, NO_VALUE, NO_PARENT, NO_FLOW,
    DisjointDenseIntSet, find_set!, union_set!, merge_a_into_b!, same_set, make_set!,
    LIFOMinPriorityQueue, pq_push!, pq_pop!,
    Depression, Outlet, dh_label_t, flat_c_idx,
    get_depression_hierarchy,
    fill_spill_merge!, depression_volume

include("dephier_oracle.jl")

const TEST_CASES_DIR = joinpath(@__DIR__, "test_cases")

const CASES = [
    (name = "case_01_trough",                ocean_level = -100.0, swl = 1.0, tol = 1e-9),
    (name = "case_02_single_depression",     ocean_level = 0.0,    swl = 1.0, tol = 1e-9),
    (name = "case_03_two_depressions",       ocean_level = 0.0,    swl = 1.5, tol = 1e-9),
    (name = "case_04_nested_metadepression", ocean_level = 0.0,    swl = 1.0, tol = 1e-9),
    (name = "case_05_perlin_100x100",        ocean_level = 0.0,    swl = 0.5, tol = 1e-6),
]

# Flip to true once `compute_julia_wtd` below is wired up to the port.
const PORT_READY = true

# Read a GeoTIFF, converting NoData cells to NaN. For input topos this gives
# the ocean border as NaN; for fsm.exe reference outputs it's a no-op since
# the binary writes 0 (not -9999) at ocean cells.
function read_tif(path::AbstractString)
    ArchGDAL.read(path) do dataset
        band = ArchGDAL.getband(dataset, 1)
        arr = ArchGDAL.read(band)
        nodata = ArchGDAL.getnodatavalue(band)
        if nodata !== nothing && !isnan(nodata)
            arr = replace(arr, nodata => NaN)
        end
        return arr
    end
end

# Single seam between this test harness and the port. Mirrors the C++
# `main.cpp` setup: label oceans, init wtd to swl on land / 0 on ocean,
# run dephier, then run FSM.
function compute_julia_wtd(topo, ocean_level, swl)
    label, flowdirs, topo_clean = prepare_label_and_flowdirs(topo, ocean_level)
    deps = get_depression_hierarchy(topo_clean, label, flowdirs)
    W, H = size(topo)
    wtd = fill(Float64(swl), W, H)
    @inbounds for i in eachindex(label)
        if label[i] == OCEAN
            wtd[i] = 0.0
        end
    end
    fill_spill_merge!(topo_clean, label, flowdirs, deps, wtd)
    return wtd
end

@testset "FillSpillMerge" begin
    @testset "constants" begin
        @testset "D8" begin
            # Each direction's inverse offsets sum to (0, 0). Catches a typo
            # in any of D8X, D8Y, or D8_INVERSE that breaks mutual consistency.
            for n in 1:8
                @test (D8X[n] + D8X[D8_INVERSE[n]], D8Y[n] + D8Y[D8_INVERSE[n]]) == (0, 0)
            end
            pairs = collect(zip(D8X, D8Y))
            @test length(unique(pairs)) == 8
            @test all(p -> p != (0, 0) && all(c -> c in (-1, 0, 1), p), pairs)
            @test length(D8X) == length(D8Y) == length(D8_INVERSE) == 8
        end
        @testset "sentinels" begin
            @test OCEAN     isa UInt32
            @test NO_DEP    isa UInt32
            @test NO_VALUE  isa UInt32
            @test NO_PARENT isa UInt32
            @test NO_FLOW   isa Int8
            # Three "max" sentinels are intentionally the same value; if a
            # future cleanup collapses them, this test fails loudly.
            @test NO_DEP === NO_VALUE === NO_PARENT
            @test NO_DEP == typemax(UInt32)
            @test OCEAN == 0
        end
    end

    @testset "fp_compare" begin
        @testset "fp_eq" begin
            @test fp_eq(1.0, 1.0)
            @test fp_eq(1.0, 1.0 + 1e-7)
            @test fp_eq(1.0, 1.0 - 1e-7)
            @test !fp_eq(1.0, 1.0 + 1e-5)
            @test !fp_eq(1.0, 1.0 - 1e-5)
            @test fp_eq(0, 0.0)
            # Strict-< boundary: exactly at the tolerance is *not* equal.
            @test !fp_eq(0.0, 1e-6)
        end
        @testset "fp_le" begin
            @test fp_le(1.0, 2.0)
            @test fp_le(1.0, 1.0)
            @test fp_le(1.0, 1.0 - 1e-7)
            @test !fp_le(2.0, 1.0)
            @test !fp_le(1.0, 1.0 - 1e-5)
        end
        @testset "fp_ge" begin
            @test fp_ge(2.0, 1.0)
            @test fp_ge(1.0, 1.0)
            @test fp_ge(1.0, 1.0 + 1e-7)
            @test !fp_ge(1.0, 2.0)
            @test !fp_ge(1.0, 1.0 + 1e-5)
        end
        @testset "symmetry" begin
            for (a, b) in [(1.0, 1.0), (1.0, 1.0 + 1e-7), (1.0, 1.0 + 1e-5), (1.0, 2.0)]
                @test fp_eq(a, b) == fp_eq(b, a)
            end
        end
        @testset "fp_eq implies fp_le and fp_ge" begin
            for (a, b) in [(1.0, 1.0), (1.0, 1.0 + 1e-7), (1.0, 1.0 - 1e-7)]
                @test fp_eq(a, b)
                @test fp_le(a, b)
                @test fp_ge(a, b)
            end
        end
        @testset "NaN" begin
            # Mirrors C++ behaviour: std::abs(NaN - x) = NaN, and NaN < 1e-6 is
            # false, so no comparison involving NaN returns true.
            @test !fp_eq(NaN, NaN)
            @test !fp_eq(NaN, 0.0)
            @test !fp_le(NaN, 0.0)
            @test !fp_ge(NaN, 0.0)
        end
    end

    @testset "bucket_fill" begin
        @testset "island (5x5)" begin
            # Outer ring at elev 0 (reached by edge fill), inner ring at elev 5
            # (the island), centre cell at elev 0 but isolated by the 5-ring.
            W, H = 5, 5
            topo = fill(5.0, W, H)
            topo[1, :] .= 0.0
            topo[W, :] .= 0.0
            topo[:, 1] .= 0.0
            topo[:, H] .= 0.0
            topo[3, 3] = 0.0

            label = fill(NO_DEP, W, H)
            bucket_fill_from_edges!(topo, label, 0.0, OCEAN)

            expected = fill(NO_DEP, W, H)
            expected[1, :] .= OCEAN
            expected[W, :] .= OCEAN
            expected[:, 1] .= OCEAN
            expected[:, H] .= OCEAN
            @test label == expected
            @test label[3, 3] == NO_DEP
        end
        @testset "empty seeds" begin
            check = zeros(5, 5)
            set = fill(NO_DEP, 5, 5)
            bucket_fill!(check, set, 0.0, OCEAN, Tuple{Int,Int}[])
            @test all(==(NO_DEP), set)
        end
        @testset "seeds on non-matching cells" begin
            check = fill(5.0, 5, 5)
            set = fill(NO_DEP, 5, 5)
            bucket_fill_from_edges!(check, set, 0.0, OCEAN)
            @test all(==(NO_DEP), set)
        end
        @testset "NaN check cells" begin
            check = fill(NaN, 5, 5)
            set = fill(NO_DEP, 5, 5)
            bucket_fill_from_edges!(check, set, NaN, OCEAN)
            @test all(==(NO_DEP), set)
        end
        @testset "exact equality (not ≤)" begin
            # Cells arbitrarily close to check_value are not labeled.
            check = fill(1e-10, 5, 5)
            set = fill(NO_DEP, 5, 5)
            bucket_fill_from_edges!(check, set, 0.0, OCEAN)
            @test all(==(NO_DEP), set)
        end
        @testset "idempotency" begin
            W, H = 5, 5
            check = fill(5.0, W, H)
            check[1, :] .= 0
            check[W, :] .= 0
            check[:, 1] .= 0
            check[:, H] .= 0
            set = fill(NO_DEP, W, H)
            bucket_fill_from_edges!(check, set, 0.0, OCEAN)
            snapshot = copy(set)
            bucket_fill_from_edges!(check, set, 0.0, OCEAN)
            @test set == snapshot
        end
        @testset "DimensionMismatch" begin
            @test_throws DimensionMismatch bucket_fill_from_edges!(zeros(5, 5), fill(NO_DEP, 4, 5), 0.0, OCEAN)
        end
        @testset "bucket_fill! with custom seed" begin
            # 3x3 patch in middle at elev 0, isolated by 5-cells on all sides.
            check = fill(5.0, 5, 5)
            for x in 2:4, y in 2:4
                check[x, y] = 0.0
            end
            set = fill(NO_DEP, 5, 5)
            bucket_fill!(check, set, 0.0, OCEAN, Tuple{Int,Int}[(3, 3)])
            expected = fill(NO_DEP, 5, 5)
            for x in 2:4, y in 2:4
                expected[x, y] = OCEAN
            end
            @test set == expected
        end
    end

    @testset "DisjointDenseIntSet" begin
        @testset "preallocated" begin
            s = DisjointDenseIntSet(5)
            for i in 0:4
                @test find_set!(s, i) == i
                @test same_set(s, i, i)
            end
        end
        @testset "union_set! by rank" begin
            s = DisjointDenseIntSet(4)
            union_set!(s, 0, 1)
            @test same_set(s, 0, 1)
            @test !same_set(s, 0, 2)
            union_set!(s, 2, 3)
            union_set!(s, 1, 3)
            for i in 0:3, j in 0:3
                @test same_set(s, i, j)
            end
        end
        @testset "merge_a_into_b! preserves parenthood" begin
            s = DisjointDenseIntSet(3)
            merge_a_into_b!(s, 0, 2)
            merge_a_into_b!(s, 1, 2)
            @test find_set!(s, 0) == 2
            @test find_set!(s, 1) == 2
            @test find_set!(s, 2) == 2
        end
        @testset "dynamic growth via make_set!" begin
            s = DisjointDenseIntSet()
            make_set!(s, 7)
            @test find_set!(s, 7) == 7
            @test_throws ArgumentError find_set!(s, 8)
        end
        @testset "merge_a_into_b! grows storage" begin
            s = DisjointDenseIntSet(2)
            merge_a_into_b!(s, 0, 5)
            @test find_set!(s, 0) == 5
            @test find_set!(s, 5) == 5
        end
    end

    @testset "LIFOMinPriorityQueue" begin
        @testset "min-key ordering" begin
            q = LIFOMinPriorityQueue{Int, Float64}()
            pq_push!(q, 3.0, 30)
            pq_push!(q, 1.0, 10)
            pq_push!(q, 2.0, 20)
            @test pq_pop!(q) == 10
            @test pq_pop!(q) == 20
            @test pq_pop!(q) == 30
            @test isempty(q)
        end
        @testset "LIFO on equal keys" begin
            # The dephier algorithm relies on this: when several cells of
            # the same elevation are queued, the most recently pushed is
            # popped first, so flat-area wavefronts stay coherent.
            q = LIFOMinPriorityQueue{Int, Float64}()
            for v in 1:5
                pq_push!(q, 1.0, v)
            end
            @test [pq_pop!(q) for _ in 1:5] == [5, 4, 3, 2, 1]
        end
        @testset "interleaved keys" begin
            q = LIFOMinPriorityQueue{Int, Float64}()
            pq_push!(q, 1.0, 1)
            pq_push!(q, 2.0, 2)
            pq_push!(q, 1.0, 3)
            pq_push!(q, 2.0, 4)
            @test pq_pop!(q) == 3   # most recent at key=1
            @test pq_pop!(q) == 1
            @test pq_pop!(q) == 4   # most recent at key=2
            @test pq_pop!(q) == 2
        end
    end

    @testset "Depression/Outlet" begin
        @testset "Depression defaults" begin
            d = Depression{Float64}()
            @test d.pit_cell == NO_VALUE
            @test d.out_cell == NO_VALUE
            @test d.parent   == NO_PARENT
            @test d.odep     == NO_VALUE
            @test d.geolink  == NO_VALUE
            @test d.pit_elev == Inf
            @test d.out_elev == Inf
            @test d.lchild   == NO_VALUE
            @test d.rchild   == NO_VALUE
            @test d.ocean_parent == false
            @test d.ocean_linked == dh_label_t[]
            @test d.dep_label == 0
            @test d.cell_count == 0
            @test d.dep_vol == 0.0
            @test d.water_vol == 0.0
            @test d.total_elevation == 0.0
        end
        @testset "Outlet swaps to depa<=depb" begin
            o1 = Outlet{Float64}(7, 3, 100, 5.0)
            @test o1.depa == 3
            @test o1.depb == 7
            o2 = Outlet{Float64}(3, 7, 100, 5.0)
            @test o2.depa == 3
            @test o2.depb == 7
        end
    end

    @testset "FSM helpers" begin
        @testset "depression_volume" begin
            # Simple sill: 3 cells at elev 1.0 dammed at elev 2.0
            # -> 3*2 - (1+1+1) = 3.
            @test depression_volume(2.0, 3, 3.0) == 3.0
            # Empty depression has zero volume regardless of sill.
            @test depression_volume(10.0, 0, 0.0) == 0.0
            # A sill at the floor of a single cell holds nothing.
            @test depression_volume(1.0, 1, 1.0) == 0.0
            # Negative-volume case (current_volume goes negative when we
            # climb past a saddle): single cell, sill below floor.
            @test depression_volume(0.0, 1, 1.0) == -1.0
        end

        @testset "fill_spill_merge! preserves ocean cells at 0" begin
            # Trivial 3x3 grid: centre is a pit at elev 1, surrounded
            # by ocean at elev 0. With swl=1, the centre receives water
            # but cannot overflow (parent is ocean-linked) so it just
            # fills to the sill (elev 0... actually it cannot fill at
            # all since dep_vol = 0). What we care about here is the
            # invariant: ocean cells stay at wtd=0.
            topo = [0.0 0.0 0.0;
                    0.0 1.0 0.0;
                    0.0 0.0 0.0]
            label, flowdirs, topo_clean = prepare_label_and_flowdirs(topo, 0.0)
            deps = get_depression_hierarchy(topo_clean, label, flowdirs)
            W, H = size(topo)
            wtd = fill(1.0, W, H)
            for i in eachindex(label)
                if label[i] == OCEAN
                    wtd[i] = 0.0
                end
            end
            fill_spill_merge!(topo_clean, label, flowdirs, deps, wtd)
            for i in eachindex(label)
                if label[i] == OCEAN
                    @test wtd[i] == 0.0
                end
            end
        end

        @testset "fill_spill_merge! is idempotent under swl=0" begin
            # With zero surface water everywhere, FSM should leave wtd
            # unchanged at zero across the board.
            topo = [0.0 0.0 0.0;
                    0.0 1.0 0.0;
                    0.0 2.0 0.0;
                    0.0 1.0 0.0;
                    0.0 0.0 0.0]
            label, flowdirs, topo_clean = prepare_label_and_flowdirs(topo, 0.0)
            deps = get_depression_hierarchy(topo_clean, label, flowdirs)
            W, H = size(topo)
            wtd = zeros(W, H)
            fill_spill_merge!(topo_clean, label, flowdirs, deps, wtd)
            @test all(==(0.0), wtd)
        end
    end

    @testset "Phase 2 oracle: $(case.name)" for case in CASES
        case_dir = joinpath(TEST_CASES_DIR, case.name)
        oracle_path = joinpath(case_dir, "expected-dh.txt")
        @test isfile(oracle_path)

        topo = read_tif(joinpath(case_dir, "input.tif"))
        oracle = read_dh_dump(oracle_path)
        @test size(topo) == (oracle.W, oracle.H)

        label, flowdirs, topo_clean = prepare_label_and_flowdirs(topo, case.ocean_level)
        deps = get_depression_hierarchy(topo_clean, label, flowdirs)

        @testset "label grid" begin
            ok, n, first = diff_grid("label", label, oracle.label)
            ok || @info "first label mismatch" first n
            @test ok
        end
        @testset "flowdirs grid" begin
            ok, n, first = diff_grid("flowdirs", flowdirs, oracle.flowdirs)
            ok || @info "first flowdirs mismatch" first n
            @test ok
        end
        @testset "depression count" begin
            @test length(deps) == length(oracle.deps)
        end
        @testset "per-depression fields" begin
            n_check = min(length(deps), length(oracle.deps))
            for i in 1:n_check
                mismatches = diff_depression(deps[i], oracle.deps[i])
                if !isempty(mismatches)
                    @info "depression mismatch" i mismatches
                end
                @test isempty(mismatches)
            end
        end
    end

    @testset "Oracle wtd: $(case.name)" for case in CASES
        case_dir = joinpath(TEST_CASES_DIR, case.name)
        @test isfile(joinpath(case_dir, "input.tif"))
        @test isfile(joinpath(case_dir, "expected-wtd.tif"))

        if PORT_READY
            topo = read_tif(joinpath(case_dir, "input.tif"))
            expected_wtd = read_tif(joinpath(case_dir, "expected-wtd.tif"))
            wtd = compute_julia_wtd(topo, case.ocean_level, case.swl)
            @test size(wtd) == size(expected_wtd)
            @test maximum(abs.(wtd .- expected_wtd)) < case.tol
        else
            @test_skip compute_julia_wtd(nothing, case.ocean_level, case.swl)
        end
    end
end
