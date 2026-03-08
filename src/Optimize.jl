module Optimize

# Write your package code here.
    export bisectv, brentqv, gaussian_quad, brentqv_parallel
    include("./bisect.jl")
    include("./brentq.jl")
    include("./gaussian_quad.jl")
end
