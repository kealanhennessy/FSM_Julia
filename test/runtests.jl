using Test
using FillSpillMerge
using ArchGDAL

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

function read_tif(path::AbstractString)
    ArchGDAL.read(path) do dataset
        return ArchGDAL.read(ArchGDAL.getband(dataset, 1))
    end
end

# Single seam between this test harness and the port. Replace the body with
# a call into the package's public FSM entry point once Phase 3 lands.
function compute_julia_wtd(topo, ocean_level, swl)
    error("Julia port not implemented yet")
end

@testset "FillSpillMerge" begin
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
