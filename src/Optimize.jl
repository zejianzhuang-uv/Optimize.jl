module Optimize

# Write your package code here.
    export bisectv, brentqv
    include("./bisect.jl")
    include("./brentq.jl")
end
