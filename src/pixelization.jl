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
"""
struct RingPixelization{R} <: AbstractRingPixelization{R}
    rings::Vector{Ring{R}}
end
# Must implement `rings`. For the rest of the interface, use the generic implementations.
rings(ringpix::RingPixelization) = ringpix.rings

# Convert any other ring map description by collecting the ring descriptions
RingPixelization(ringpix::AbstractRingPixelization) = RingPixelization(collect(rings(ringpix)))
