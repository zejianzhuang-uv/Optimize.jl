module Optimize

# Write your package code here.
    export bisectv, brentqv, gaussian_quad
    include("./bisect.jl")
    include("./brentq.jl")
    include("./gaussian_quad.jl")
end
