using Optimize
using Test

@testset "Optimize.jl" begin

    derivative(x -> exp(x), 2.) |> println
end
