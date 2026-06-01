using Optimizer
using Test

@testset "Optimizer.jl" begin

    derivative(x -> exp(x), 2.) |> println
end
