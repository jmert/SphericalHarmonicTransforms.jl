using AssociatedLegendrePolynomials
using FFTW
@static if VERSION < v"1.6.0-DEV.1591"
    using Compat: cispi # Compat@v3.25
end

# Brute-force synthesis of the real-space points at (θ, ϕ) given a set of alms.
# The idea is that there are no optimizations, so it is much easier to verify the
# correctness of the function at the expense of a huge hit to performance.
function synthesize_reference(alms::AbstractMatrix{C}, θ, ϕ) where {C<:Complex}
    axes(θ) == axes(ϕ) ||
        throw(DimensionMismatch("`θ` and `ϕ` must have same axes"))
    lmax, mmax = size(alms) .- 1

    R = real(C)
    Λ = Matrix{R}(undef, size(alms)...)
    Φ = Vector{C}(undef, mmax + 1)
    syn = zeros(R, axes(θ)...)

    for I in eachindex(syn)
        λlm!(Λ, lmax, mmax, cos(θ[I]))       # λ_ℓ^m(cos θ) factors
        @. Φ = cispi.((0:mmax) .* (ϕ[I]/π))  # e^{imϕ} factors
                                             #   Using cispi(ϕ/π) rather than cis(ϕ) gives
                                             #   slightly more accurate results for test
                                             #   healpix rings.

        acc = zero(R)
        # Σ_{ℓ = 0}^{ℓmax}
        for ℓ in 0:lmax
            # Σ_{m = 0}^{mmax}
            acc += real(alms[ℓ+1,1] * Λ[ℓ+1,1])
            for m in 1:min(ℓ, mmax)
                # Assuming alms were sourced from real field, alms are constrained such
                # that
                #    a_{ℓ(-m)} Y_ℓ^(-m) + a_{ℓ(+m)} Y_ℓ^(+m)
                #    == 2 * Re[ a_{ℓ(+m)} Y_ℓ^(+m) ]
                acc += 2 * real(alms[ℓ+1,m+1] * Λ[ℓ+1,m+1] * Φ[m+1])
            end
        end
        syn[I] = acc
    end
    return syn
end

# Similar to above, a brute-force analysis of the given map (with corresponding coordinates
# (θ, ϕ) for each pixel into the harmonic coefficients up to order lmax and degree mmax.
#
# Note that there is no information provided on the area (in terms of the integral's
# surface area element) of pixels, so that normalization factor is not applied during
# analysis and must be applied by the caller using information it will have.
#
# For instance, an ECP grid needs an extra normalization factor of
#
#   alms .*= 2π/N_ϕ * π/N_θ
#        .*= 2π²/N_pix
#
# Without quadrature weights, an iterative approach will be necessary to achieve any
# kind of relative accuracy:
#
#   lmax = size(alms, 1) - 1
#   map = synthesize_ecp(alms, Nθ, Nϕ)
#   θ, ϕ = make_grid(Nθ, Nϕ)
#
#   alms′ = (2π^2 / prod(size(map))) * analyze_reference(map, θ, ϕ, lmax)
#   map′ = synthesize_ecp(alms′, Nθ, Nϕ)
#   δalms = (2π^2 / prod(size(map))) * analyze_reference(map - map′, θ, ϕ, lmax)
#   alms += δalms  # iterative refinement
function analyze_reference(map, θ, ϕ, lmax::Integer, mmax::Integer = lmax)
    axes(map) == axes(θ) == axes(ϕ) ||
        throw(DimensionMismatch("`map`, `θ`, and `ϕ` must all have the same axes"))

    R = eltype(map)
    C = complex(R)
    Λ = Matrix{R}(undef, lmax + 1, mmax + 1)
    Φ = Vector{C}(undef, mmax + 1)
    alms = zeros(C, lmax + 1, mmax + 1)

    for I in eachindex(map)
        sθ, cθ = sincos(θ[I])
        λlm!(Λ, lmax, mmax, cθ)             # λ_ℓ^m(cos θ) factors
        @. Φ = cispi((0:mmax) * (-ϕ[I]/π))  # e^{-imϕ} factors
                                            #   Using cispi(ϕ/π) rather than cis(ϕ) gives
                                            #   slightly more accurate results for test
                                            #   healpix rings.

        # Σ_{ℓ = 0}^{ℓmax}
        for ℓ in 0:lmax
            # Σ_{m = 0}^{mmax}
            alms[ℓ+1,1] += map[I] * Λ[ℓ+1,1] * sθ
            for m in 1:min(ℓ, mmax)
                alms[ℓ+1,m+1] += map[I] * Λ[ℓ+1,m+1] * Φ[m+1] * sθ
            end
        end
    end
    return alms
end
