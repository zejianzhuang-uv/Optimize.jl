

# function derivative(f::Function, x0::Union{Float64, ComplexF64}; eps_start=1e-4, tol=1e-6, max_iter=20)
#     """
#     Numerical derivative with convergence check.
#     Returns derivative value. Prints warning if not converged.
#     """
#     eps = eps_start
#     deri = (f(x0 + eps) - f(x0)) / eps
    
#     for iter in 1:max_iter
#         deri_new = (f(x0 + eps) - f(x0)) / eps
#         error = abs(deri_new - deri)
        
#         if error < tol
#             return deri_new
#         end
        
#         deri = deri_new
#         eps /= 10
#     end
    
#     @warn "Derivative did not converge after $max_iter iterations. Last value: $deri"
#     return deri
# end

function derivative(f::Function, x0::Union{Float64,ComplexF64}; tol=1e-8, max_iter=50)
    if isnan(x0)
        return NaN64
    end
    scale = 1e12
    h = eps() * scale
    deri = (f(x0 + h) - f(x0)) / h
    iter = 1
    while iter < max_iter
        scale /= 10
        h = eps(typeof(x0)) * scale
        deri_new = (f(x0 + h) - f(x0)) / h
        error = abs(deri_new - deri)
        deri = deri_new
        iter += 1

        if error < tol
            return deri_new
        end
    end

    @warn "Not converged. Last h: $h"
    return deri
end