

function derivative(f::Function, x0::Real; eps_start=1e-4, tol=1e-6, max_iter=20)
    """
    Numerical derivative with convergence check.
    Returns derivative value. Prints warning if not converged.
    """
    eps = eps_start
    deri = (f(x0 + eps) - f(x0)) / eps
    
    for iter in 1:max_iter
        deri_new = (f(x0 + eps) - f(x0)) / eps
        error = abs(deri_new - deri)
        
        if error < tol
            return deri_new
        end
        
        deri = deri_new
        eps /= 10
    end
    
    @warn "Derivative did not converge after $max_iter iterations. Last value: $deri"
    return deri
end