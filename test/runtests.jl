using Test
using Logging
using Statistics
using Printf
using JLD2
using CUDA

using Oceananigans
using Oceananigans.Architectures
using JULES

Logging.global_logger(OceananigansLogger())

Archs = CUDA.has_cuda() ? [CPU] : [CPU]

@testset "JULES" begin
    include("test_models.jl")
    include("test_lazy_fields.jl")
    include("test_time_stepping.jl")
    include("test_regression.jl")
end
