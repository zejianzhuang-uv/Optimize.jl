module Optimize

# Write your package code here.
    export bisectv, brentqv, gaussian_quad, brentqv_parallel, derivative
    include("./bisect.jl")
    include("./brentq.jl")
    include("./gaussian_quad.jl")
    include("./derivative.jl")
end
