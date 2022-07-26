"""
    Ring{R<:Real}

Parameters required to describe an isolatitude ring on the sphere compatible with spherical
harmonic transforms.

### Fields
- `offset::Tuple{Int, Int}`:
  Offsets (1-indexed) to the first pixel of the isolatitude ring(s) within a map array.
  The first offset in the tuple corresponds to the colatitude ``θ`` given by `.cosθ`, while
  the second offset is for the ring mirrored over the equator (colatitude ``π - θ``).
  Either offset may be 0 to indicate the ring is not present (such as for the equator, but
  also partial sphere coverage such as transforms on bands).
- `stride::Int`: Stride between pixels within the ring.
- `nϕ::Int`: The number of equispaced pixels in the isolatitude ring.
- `cosθ::R`: Cosine of the colatitude of the ring, i.e. ``\\cos(θ)``.
- `ϕ₀_π::R`: Azimuth of the first pixel in the isolatitude ring ``ϕ₀`` divided by ``π``.
- `ΔΩ::R`: Solid angle area of pixels in the ring, roughly ``\\sin(θ) × Δθ × Δϕ``.
"""
struct Ring{R}
    offset::Tuple{Int,Int}
    stride::Int
    nϕ::Int
    cosθ::R
    ϕ₀_π::R
    ΔΩ::R
end

"""
    abstract type AbstractRingPixelization{R} end

Abstract supertype of all ring map descriptor types.
"""
abstract type AbstractRingPixelization{R} end

_pixeltype(::AbstractRingPixelization{R}) where {R} = R

function Base.show(io::IO, ::MIME"text/plain", ringpix::AbstractRingPixelization)
    if get(io, :compact, false)
        show(io, ringpix)
    else
        props = _summary(ringpix)
        print(io, props.nring, "-ring ", summary(ringpix), " with ", props.npix, " pixels")
    end
end

"""
    rr = rings(ringpix::AbstractRingPixelization)

A vector (or iterable) `rr` of [`Ring`](@ref)s described by `ringpix`.
"""
function rings end

function _summary(ringpix::AbstractRingPixelization)
    nring, npix, nϕmax = 0, 0, typemin(Int)
    for rinfo in rings(ringpix)
        o₁, o₂ = rinfo.offset
        mult = (o₁ != 0) + (o₂ != 0)
        nring += mult
        npix += mult * rinfo.nϕ
        nϕmax = mult != 0 ? max(nϕmax, rinfo.nϕ) : nϕmax
    end
    return (; nring = nring, npix = npix, nϕmax = nϕmax)
end

"""
    nr = nring(ringpix::AbstractRingPixelization)

Number of isolatitude rings `nr` described by `ringpix`.
"""
function nring(ringpix::AbstractRingPixelization)
    return _summary(ringpix)[:nring]
end

"""
    np = npix(ringpix::AbstractRingPixelization)

Number of pixels `np` covered by all rings described by `ringpix`.
"""
function npix(ringpix::AbstractRingPixelization)
    return _summary(ringpix)[:npix]
end

"""
    n = nϕmax(ringpix::AbstractRingPixelization)

Number of pixels `n` in the longest isolatitude ring described by `ringpix`.
"""
function nϕmax(ringpix::AbstractRingPixelization)
    return _summary(ringpix)[:nϕmax]
end

"""
    mapbuf = buffer(ringpix::AbstractRingPixelization, ::Type{T})

Allocates an array `mapbuf` with element type `T` appropriate for containing a map in the
pixelization scheme described by `ringpix`.
"""
buffer(ringpix::AbstractRingPixelization, ::Type{T}) where {T} = Vector{T}(undef, npix(ringpix))
buffer(ringpix::AbstractRingPixelization) = buffer(ringpix, _pixeltype(ringpix))


"""
    RingPixelization{R} <: AbstractRingPixelization{R}

A pixelization of the unit sphere described by an arbitrary vector of [`Ring`](@ref)s.

# Examples
```jldoctest; setup = :(import SphericalHarmonicTransforms: ECPPixelization, RingPixelization)
julia> ringpix = RingPixelization(ECPPixelization(5, 10))
5-ring RingPixelization{Float64} with 50 pixels
```
"""
struct RingPixelization{R} <: AbstractRingPixelization{R}
    rings::Vector{Ring{R}}
end
# Must implement `rings`. For the rest of the interface, use the generic implementations.
rings(ringpix::RingPixelization) = ringpix.rings

# Convert any other ring map description by collecting the ring descriptions
RingPixelization(ringpix::AbstractRingPixelization) = RingPixelization(collect(rings(ringpix)))


@doc raw"""
    ECPPixelization{R} <: AbstractRingPixelization{R}

An Equidistant Cylindrical Projection pixelization in column-major order, with dimensions
`nθ × nϕ` in the colatitude/azimuth directions, respectively.

See also [`RingPixelization`](@ref)

# Examples
```jldoctest; setup = :(import SphericalHarmonicTransforms: ECPPixelization)
julia> ECPPixelization(250, 500)
250×500 ECPPixelization{Float64}

julia> ECPPixelization{Float32}(250, 500)
250×500 ECPPixelization{Float32}
```

# Extended help

The described pixelization covers the entire sphere ``[0, π] × [0, 2π]`` with
uniformly-sized pixels in _coordinate space_ of size ``Δθ = π/n_θ`` by ``Δϕ = 2π/n_ϕ``.
The pixel _centers_ take on coordinates
```math
    (θ_j, ϕ_k) = \left( π \frac{2j + 1}{2n_θ}, 2π \frac{2k + 1}{2n_ϕ} \right)
```
for integers ``j ∈ \{0, …, n_θ - 1\}`` and ``k ∈ \{0, …, n_ϕ - 1\}``.
On the sphere, the pixels cover a solid angle ``ΔΩ_{jk} ≈ \sin(θ_j) Δθ Δϕ``, i.e. the
pixels' physical size approaches zero towards the poles.

See also
[equidistant cylindrical projection](https://en.wikipedia.org/wiki/Equirectangular_projection)
"""
struct ECPPixelization{R} <: AbstractRingPixelization{R}
    nθ::Int
    nϕ::Int
end
ECPPixelization(nθ, nϕ) = ECPPixelization{Float64}(Int(nθ), Int(nϕ))

function Base.show(io::IO, ::MIME"text/plain", ecppix::ECPPixelization)
    if get(io, :compact, false)
        show(io, ecppix)
    else
        print(io, ecppix.nθ, "×", ecppix.nϕ, " ", summary(ecppix))
    end
end

function rings(ecppix::ECPPixelization)
    T = _pixeltype(ecppix)
    nθ, nϕ = ecppix.nθ, ecppix.nϕ
    # rings are symmetric over the equator, so only require half
    nθh = (nθ + 1) ÷ 2
    # ϕ offset and ΔθΔϕ are same for all rings
    ϕ_π = one(T) / nϕ
    ΔθΔϕ = 2T(π)^2 / (nθ * nϕ)
    function ecpring(ii)
        o₁ = ii
        o₂ = nθ - ii + 1
        offs = (o₁, isodd(nθ) && ii == nθh ? 0 : o₂)
        sθ, cθ = sincospi(T(2ii - 1) / 2nθ)
        return Ring{T}(offs, nθ, nϕ, cθ, ϕ_π, sθ * ΔθΔϕ)
    end
    return Broadcast.instantiate(Broadcast.broadcasted(ecpring, 1:nθh))
end

# overload to return square map
buffer(ecppix::ECPPixelization, ::Type{T}) where {T} = Matrix{T}(undef, ecppix.nθ, ecppix.nϕ)
# optimized overloads
nring(ecppix::ECPPixelization) = ecppix.nθ
npix(ecppix::ECPPixelization) = ecppix.nθ * ecppix.nϕ
nϕmax(ecppix::ECPPixelization) = ecppix.nϕ
