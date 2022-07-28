import Compat.@compat  # Compat@v3.21 for @compat import Mod as NewName
@compat import AssociatedLegendrePolynomials as Legendre

using .Legendre: λlm!, unsafe_legendre!
using Base: @propagate_inbounds
using FFTW
using LinearAlgebra: mul!

@static if VERSION < v"1.6.0-DEV.1591"
    using Compat: sincospi # Compat@v3.23
    using Compat: cispi # Compat@v3.25
end


"""
    alias_index(len::Int, m::Int) -> (i::Int, isconj::Bool, isnyq::Bool)

Aliases the given 0-indexed `m` for a periodic FFT frequency axis of length `len` --- with
an additional aliasing to below the Nyquist frequency approprate for use in a real-only
[`irfft()`] transform --- to the new index `i`.

To handle symmetries of the FFT, the flags `isconj` and `isnyq` indicate if the Fourier
coefficients at `i` should be conjugated or real-doubled (i.e. `c = complex(2real(c))`),
respectively.
"""
@inline function alias_index(len::Int, i::Int)
    isconj, isnyq = false, false
    nyq = max(1, len ÷ 2)
    i < nyq && return (i, isconj, isnyq)
    i = mod(i, len)
    if i > nyq
        i = len - i
        isconj = true
    elseif i == 0 || (iseven(len) && i == nyq)
        isnyq = true
    end
    return (i, isconj, isnyq)
end

"""
    alias_coeffs(coeffs, isconj::Bool, isnyq::Bool) -> coeffs

Helper for [`alias_index`](@ref) which applies the conjugation or pure-real symmetry
transformations to the Fourier coefficients `coeffs` based on the two flags `isconj` and
`isnyq`.
"""
@inline function alias_coeffs(coeffs, isconj::Bool, isnyq::Bool)
    return ifelse(isnyq, complex.(2 .* real.(coeffs)),
                  ifelse(isconj, conj.(coeffs), coeffs))
end

function _promoteeltype(ringpix::AbstractRingPixelization, arrs...)
    P = _pixeltype(ringpix)
    Ts = real.(eltype.(arrs))
    R = promote_type(P, Ts...)
    C = complex(R)
    return R, C
end

"""
    mapbuf = synthesize(ringpix::AbstractRingPixelization, alms::AbstractMatrix{<:Complex})

Performs the spherical harmonic synthesis transform of the harmonic coefficients `alms`
for a real-valued map `mapbuf`, pixelized as described by `ringpix`.

See also [`synthesize!`](@ref), [`analyze`](@ref), [`analyze!`](@ref)
"""
function synthesize(ringpix::AbstractRingPixelization, alms::AbstractMatrix{<:Complex})
    R, C = _promoteeltype(ringpix, alms)
    return synthesize!(buffer(ringpix, R), ringpix, alms)
end

@doc raw"""
    synthesize!(mapbuf::AbstractArray{<:Real}, ringpix::AbstractRingPixelization, alms::AbstractMatrix{<:Complex})

Performs the spherical harmonic synthesis transform of the harmonic coefficients `alms`
for a real-valued map `mapbuf`, pixelized as described by `ringpix`.

See also [`synthesize`](@ref), [`analyze`](@ref), [`analyze!`](@ref)

# Extended help

The synthesis transform is implemented such that a real-valued function ``f(θ, ϕ)``
is generated from complex-valued spherical harmonic coefficients ``a_{ℓm}`` as
```math
    f(θ, ϕ) = \sum_{ℓ=0}^{ℓ_\mathrm{max}} \sum_{m=0}^{ℓ} a_{ℓm} Y_{ℓm}(θ, ϕ)
```
where the ``a_{ℓm}``s matrix contains only the non-negative orders (``m ≥ 0``) since the
negative orders are not independent when ``f(θ, ϕ)`` is real.
"""
function synthesize!(mapbuf::AbstractArray{<:Real}, ringpix::AbstractRingPixelization, alms::AbstractMatrix{<:Complex})
    R, C = _promoteeltype(ringpix, mapbuf, alms)
    lmax, mmax = size(alms) .- 1

    nϕ_max = nϕmax(ringpix)
    rv = Vector{R}(undef, nϕ_max)
    f₁ = zeros(C, (nϕ_max ÷ 2) + 1)
    f₂ = zeros(C, (nϕ_max ÷ 2) + 1)

    Λ = zeros(R, lmax + 1, mmax + 1)
    Λw = Legendre.Work(λlm!, Λ, Legendre.Scalar(zero(R)))

    for rpix in rings(ringpix)
        unsafe_legendre!(Λw, Λ, lmax, mmax, rpix.cosθ)

        o₁, o₂ = rpix.offset
        nϕ = rpix.nϕ
        nϕh = nϕ ÷ 2 + 1
        f₁′ = @view f₁[1:nϕh]
        f₂′ = @view f₂[1:nϕh]
        F = plan_brfft(f₁′, nϕ)

        for m in 0:mmax
            a₁, a₂ = zero(C), zero(C)
            for ℓ in m:lmax
                term = alms[ℓ+1,m+1] * Λ[ℓ+1,m+1]
                a₁ += term
                a₂ += isodd(ℓ + m) ? -term : term
            end
            i, isconj, isnyq = alias_index(rpix.nϕ, m)
            a₁, a₂ = (a₁, a₂) .* cispi(m * rpix.ϕ₀_π)
            a₁, a₂ = alias_coeffs((a₁, a₂), isconj, isnyq)
            if o₁ != 0
                f₁′[i+1] += a₁
            end
            if o₂ != 0
                f₂′[i+1] += a₂
            end
        end

        to_nϕ = Base.OneTo(nϕ)
        @inline _strided(o) = range(o, step = rpix.stride, length = nϕ)

        rv′ = @view rv[1:nϕ]
        if o₁ != 0
            mul!(rv′, F, f₁′)
            copyto!(mapbuf, _strided(o₁), rv, to_nϕ)
            fill!(f₁′, zero(C))
        end
        if o₂ != 0
            mul!(rv′, F, f₂′)
            copyto!(mapbuf, _strided(o₂), rv, to_nϕ)
            fill!(f₂′, zero(C))
        end
    end
    return mapbuf
end

"""
    alms = analyze(ringpix::AbstractRingPixelization, mapbuf::AbstractArray{<:Real}, lmax::Integer, mmax::Integer = lmax)

Performs the spherical harmonic analysis transform of the real-valued map `mapbuf`,
pixelized as described by `ringpix`, to the corresponding complex-valued spherical
harmonic coefficients `alms`.

See also [`analyze!`](@ref), [`synthesize`](@ref), [`synthesize!`](@ref)
"""
function analyze(ringpix::AbstractRingPixelization, mapbuf::AbstractArray{<:Real}, lmax::Integer, mmax::Integer = lmax)
    R, C = _promoteeltype(ringpix, mapbuf)
    alms = zeros(C, lmax + 1, mmax + 1)
    return analyze!(alms, ringpix, mapbuf)
end

@doc raw"""
    analyze!(alms::AbstractMatrix{<:Complex}, ringpix::AbstractRingPixelization, mapbuf::AbstractAray{<:Real})

Performs the discrete spherical harmonic analysis transform of the real-valued map `mapbuf`,
pixelized as described by `ringpix`, to the corresponding complex-valued spherical
harmonic coefficients `alms`.

See also [`analyze`](@ref), [`synthesize`](@ref), [`synthesize!`](@ref)

# Extended help

The analysis transform is implemented such that a real-valued function ``f(θ, ϕ)``
is analyzed for its corresponding complex-valued spherical harmonic coefficients ``a_{ℓm}``
as
```math
    a_{ℓm} = \sum_{ℓ=0}^{ℓ_\mathrm{max}} \sum_{m=0}^ℓ f(θ, ϕ) \overline{Y_{ℓm}}(θ, ϕ) ΔΩ
```
where the overline indicates complex conjugation, ``ΔΩ`` is the integration element
(pixel area) on the sphere, and the ``a_{ℓm}``s matrix contains only the non-negative orders
(``m ≥ 0``) since the negative orders are not independent when ``f(θ, ϕ)`` is real.
"""
function analyze!(alms::AbstractMatrix{<:Complex}, ringpix::AbstractRingPixelization, mapbuf::AbstractArray{<:Real})
    R, C = _promoteeltype(ringpix, alms, mapbuf)
    lmax, mmax = size(alms) .- 1

    nϕ_max = nϕmax(ringpix)
    rv = Vector{R}(undef, nϕ_max)
    f₁ = Vector{C}(undef, (nϕ_max ÷ 2) + 1)
    f₂ = Vector{C}(undef, (nϕ_max ÷ 2) + 1)

    Λ = zeros(R, lmax + 1, mmax + 1)
    Λw = Legendre.Work(λlm!, Λ, Legendre.Scalar(zero(R)))

    for rpix in rings(ringpix)
        unsafe_legendre!(Λw, Λ, lmax, mmax, rpix.cosθ)

        o₁, o₂ = rpix.offset
        nϕ = rpix.nϕ
        nϕh = nϕ ÷ 2 + 1
        rv′ = @view rv[1:nϕ]
        F = plan_rfft(rv′, 1)

        to_nϕ = Base.OneTo(nϕ)
        @inline _strided(o) = range(o, step = rpix.stride, length = nϕ)

        f₁′ = @view f₁[1:nϕh]
        f₂′ = @view f₂[1:nϕh]
        if o₁ != 0
            copyto!(rv, to_nϕ, mapbuf, _strided(o₁))
            mul!(f₁′, F, rv′)
        end
        if o₂ != 0
            copyto!(rv, to_nϕ, mapbuf, _strided(o₂))
            mul!(f₂′, F, rv′)
        end

        for m in 0:mmax
            i, isconj, isnyq = alias_index(rpix.nϕ, m)
            a₁ = o₁ != 0 ? f₁′[i+1] : zero(C)
            a₂ = o₂ != 0 ? f₂′[i+1] : zero(C)
            a₁, a₂ = alias_coeffs((a₁, a₂), isconj, isnyq)
            a₁, a₂ = (a₁, a₂) .* (rpix.ΔΩ * cispi(m * -rpix.ϕ₀_π))
            for ℓ in m:lmax
                c = isodd(ℓ+m) ? a₁ - a₂ : a₁ + a₂
                alms[ℓ+1,m+1] += c * Λ[ℓ+1,m+1]
            end
        end
    end
    return alms
end
