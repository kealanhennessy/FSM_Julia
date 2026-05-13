using Test
using FillSpillMerge
using ArchGDAL

import FillSpillMerge:
    fp_eq, fp_le, fp_ge,
    bucket_fill!, bucket_fill_from_edges!,
    D8X, D8Y, D8_INVERSE,
    OCEAN, NO_DEP, NO_VALUE, NO_PARENT, NO_FLOW

const TEST_CASES_DIR = joinpath(@__DIR__, "test_cases")

const CASES = [
    (name = "case_01_trough",                ocean_level = -100.0, swl = 1.0, tol = 1e-9),
    (name = "case_02_single_depression",     ocean_level = 0.0,    swl = 1.0, tol = 1e-9),
    (name = "case_03_two_depressions",       ocean_level = 0.0,    swl = 1.5, tol = 1e-9),
    (name = "case_04_nested_metadepression", ocean_level = 0.0,    swl = 1.0, tol = 1e-9),
    (name = "case_05_perlin_100x100",        ocean_level = 0.0,    swl = 0.5, tol = 1e-6),
]

# Flip to true once `compute_julia_wtd` below is wired up to the port.
const PORT_READY = false

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

# Single seam between this test harness and the port. Replace the body with
# a call into the package's public FSM entry point once Phase 3 lands.
function compute_julia_wtd(topo, ocean_level, swl)
    error("Julia port not implemented yet")
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

    @testset "Oracle: $(case.name)" for case in CASES
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
