

function chi_square(yth::AbstractVector{Float64}, yexp::AbstractVector{Float64}, inv_σ::AbstractVector{Float64})
    dy = @. (yth - yexp) * inv_σ # broadcast
    chi2 = sum(x -> x^2, dy)
    return chi2
end

function chi_square(yth::AbstractVector{Float64}, yexp::AbstractVector{Float64}, inv_σ::AbstractMatrix{Float64})
    dy = yth - yexp 
    chi2 = dy' * inv_σ * dy # broadcast
    return chi2
end

function weight_chi_square(yth::AbstractVector{Float64}, yexp::AbstractVector{Float64}, inv_σ::AbstractVector{Float64}, nk::Int64)
    chi2 = chi_square(yth, yexp, inv_σ)
    return chi2 / nk
end

function weight_chi_square(yth::AbstractVector{Float64}, yexp::AbstractVector{Float64}, inv_σ::AbstractMatrix{Float64}, nk::Int64)
    chi2 = chi_square(yth, yexp, inv_σ)
    return chi2 / nk
end

"""
Weighted-reduced χ square
- Nsp: The number of the scattering processes: K⁻p → {π⁰Λ, ...} and π⁻p → {ηn, ...}
- nₖ: The number of the observables in the scattering process k
- Ntot: Total number of observables in Nsp scattering processes
"""
function weight_redχ_square(all_weight_chi2::Float64, Ntot::Int64, nfit::Int64, Nsp::Int64)
    chi2 = (Ntot / (Ntot - nfit)) * (1/Nsp) * all_weight_chi2
    return chi2
end