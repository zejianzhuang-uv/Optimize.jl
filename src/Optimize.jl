module Optimize

# Write your package code here.
    export bisectv, brentqv, gaussian_quad, brentqv_parallel, derivative, STD, df_STD
    include("./bisect.jl")
    include("./brentq.jl")
    include("./gaussian_quad.jl")
    include("./derivative.jl")
    include("STD.jl")
end
