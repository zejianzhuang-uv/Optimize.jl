module Optimizer

# Write your package code here.
include("./bisect.jl")
include("./brentq.jl")
include("./gaussian_quad.jl")
include("./derivative.jl")
include("STD.jl")
include("weight_chi_square.jl")




export bisectv, brentqv, gaussian_quad, derivative, STD, df_STD
# chi2
export chi_square, weight_chi_square, weight_redχsq
end
