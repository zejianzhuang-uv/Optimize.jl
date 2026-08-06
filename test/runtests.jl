using Optimizer
using Test

@testset "Optimizer.jl" begin

    @time display(brentqv(sin, -4π:0.1:4π, 3) )
end
