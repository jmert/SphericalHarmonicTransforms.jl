# Developer Documentation

```@meta
CurrentModule = SphericalHarmonicTransforms
DocTestFilters = Regex[
        r"Ptr{0x[0-9a-f]+}",
        r"[0-9\.]+ seconds( \(.*\))?",
        ]
```

```@contents
Pages = ["devdocs.md"]
Depth = 2
```

## [Adding New Pixelizations](@id pixelnew)

### [Necessary Properties of Ring-based Pixelizations](@id pixelproperties)

As the name "ring-based pixelizations" suggests, the first obvious requirement for supported
pixelization schemes is that the format be describable in terms of isolatitude rings.
This restriction greatly reduces the number of times the associated Legendre polynomials
must be calculated — from (in principle) every pixel if they all have a unique latitude to
only once per isolatitude ring.
Therefore, the first per-ring parameter we must be told is the colatitude ``\theta`` of each
ring.

The second requirement is that each isolatitude ring's pixels meet the necessary conditions
for the Fast Fourier Transform (FFT) to be used — namely, the pixels in azimuth are
uniformly spaced and span (exactly) a full ``2\pi`` radians.
These conditions are met when the azimuth coordinates are at locations
``\phi_k = \phi_0 + 2\pi k / N_\phi``
(for integer index ``k = [0, 1, \ldots, N_\phi - 1]``), and therefore two additional
parameters must be given:
the number of pixels within the ring ``N_\phi`` and the azimuth offset of the first pixel in
the ring ``\phi_0``.
For numerical efficiency and accuracy reasons, we will instead use the azimuth offset
divided by pi: ``\phi_0/\pi``.

These three parameters (``\theta``, ``N_\phi``, ``\phi_0``) are sufficient for implementing
the synthesis transform of a single isolatitude ring.
For the analysis transform, though, we also require the integration
[volume element](https://en.wikipedia.org/wiki/Volume_element)
``d\Omega \equiv \sin\theta \, d\theta \, d\phi``, or rather its finite approximation
``\Delta\Omega \approx \sin\theta \, \Delta\theta \, \Delta\phi``.
Taking the volume element ``\Delta\Omega`` as a required parameter, the quantity
``\sin\theta`` is never calculated explicitly;
instead only ``\cos\theta`` is used as the argument of the associated Legendre polynomials,
so we can actually also revise the first parameter to instead be the cosine of the
colatitude angle ``\cos\theta``.

Together, this gives us 4 parameters which are required to sufficiently describe the
geometry of an isolatitude ring so that we can perform spherical harmonic transforms on it:

1. The cosine of the ring's colatitude angle: ``\cos\theta``
2. The azimuth angle, divided by ``\pi``, of the first pixel in the ring: ``\phi_0/\pi``
3. The number of pixels equally spaced in azimuth within the ring: ``N_\phi``
4. The surface area of each pixel within the ring: ``\Delta\Omega``

What remains is to describe how multiple rings are assembled into a map that covers the
sphere, and here again the ECP grid proves a useful case study.
The ECP grid maps well to a matrix (2D array, where colatitude increases across rows and
azimuth across columns), but there is an inherent ambiguity in whether the ordering of the
matrix in memory is
[row-major or column-major](https://en.wikipedia.org/wiki/Row-_and_column-major_order).
If we choose one — say column-major order, consistent with
[Julia's matrix representation](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-column-major)
— then the contents of a given isolatitude ring are non-contiguous in memory.
One option would be to require that the transpose of the ECP map matrix be used (so that a
ring is a contiguous vector of elements), but another more flexible option is to instead
make the [stride](https://en.wikipedia.org/wiki/Stride_of_an_array) ``s`` of the array an
additional ring description parameter.
When paired with offset of the first pixel in the ring $o$, we have sufficient information
to completely describe an ECP grid in either row- or column-major order.

As one final optimization, though, recall that if the pixel grid is symmetric over the
equator, then rings are best considered in _pairs_, since the associated Legendre
polynomials are the same up to a negative sign for colatitude angles
``\theta$ and $\pi - \theta``.
If we replace the single parameter ``o`` with a tuple of parameters ``(o_1, o_2)``, we can
inform our implementation explicitly whether the ring-pair optimization applies (such as for
a north-south pair of rings) or not (like for the singular equatorial row) by encoding
either two or one valid indices in the tuple, respectively.

Finally, we add these three additional parameters to the above list:

5. Offset to the first pixel in the northern- and southern-hemisphere rings, respectively
   (if valid): ``(o_1, o_2)``
6. The stride across the map array between successive pixels within the ring: ``s``


### [Pixelization Interface](@id pixelinterface)

**Required**

| Interfaces to extend/implement     | Brief description                                                                                |
|:---------------------------------- |:------------------------------------------------------------------------------------------------ |
| [`AbstractRingPixelization`](@ref) | Supertype of ring-based map pixelizations                                                        |
| [`rings()`](@ref)                  | Returns an indexable container of [`Ring`](@ref) descriptions for each ring in the described map |

**Optional Overrides**

The following functions have a generic implementation which deduces the necessary values
by inspecting every ring (returned by [`rings`](@ref)), but a particular pixelization is
encouraged to specialize these if they can be calculated more efficiently using structural
knowledge of the pixelization.

| Interfaces to extend/implement | Brief description                                                                     |
|:------------------------------ |:------------------------------------------------------------------------------------- |
| [`buffer()`](@ref)             | Returns an array compatible with storing a map described by a given ring pixelization |
| [`nring()`](@ref)              | The number of isolatitude rings in the pixelization                                   |
| [`npix()`](@ref)               | The total number of pixels in the pixelization                                        |
| [`nϕmax()`](@ref)              | The number of pixels in the longest isolatitude ring in the pixelization              |

### [Example Implementation](@id pixelexample)

The built-in ECP pixelization does not support exact quadrature on the sphere, but it is
included because the axes are uniformly sampled (so the maps are easily represented as
matrices without any special handling) and the simplicity of its implementation.
For more than trivial or example cases, though, it is very likely that another pixelization
scheme will be required.

A well-known improvement is to modify the ECP grid to permit non-uniformly sampled points in
``\theta``, with the node locations (and weights, see below) corresponding to the
[Gauss-Legendre nodes](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature).
In this section, we will demonstrate implementing support for spherical harmonic transforms
on the the Gauss-Legendre (GL) pixelization.
Computation of the Gauss-Legendre nodes and weights is provided by
[`FastGaussQuadrature.jl`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl).

The first step is to import the [interface](@ref pixelinterface) types.
The new `GLPixelization` type must be a subclass of the `AbstractRingPixelization`, where
the type parameter `R` is the nominal element type used during the transform.

```jldoctest examplepix
julia> import SphericalHarmonicTransforms: AbstractRingPixelization, Ring,
               rings, buffer, nring, npix, nϕmax

julia> struct GLPixelization{R} <: AbstractRingPixelization{R}
           nθ::Int
           nϕ::Int
       end

julia> GLPixelization(nθ, nϕ) = GLPixelization{Float64}(nθ, nϕ)
GLPixelization
```

Then, the only other **required** step is to implement a specialization of the
[`rings`](@ref) method to return a vector of [`Ring`](@ref) structures that describe the
[properties of the ring](@ref pixelproperties).

```jldoctest examplepix
julia> using FastGaussQuadrature: gausslegendre

julia> function rings(glpix::GLPixelization{R}) where {R}
           nθ, nϕ = glpix.nθ, glpix.nϕ
           nθh = (nθ + 1) ÷ 2  # number of rings in N. hemisphere + equator (if nθ is odd)
           ϕ_π = one(R) / nϕ   # constant pixel offset for all rings
           Δϕ = 2R(π) / nϕ     # constant azimuthal pixel size for all rings
           stride = nθ         # constant stride, column-major ordering
           cθ, sθΔθ = gausslegendre(nθ)  # z = cos(θ) and corresponding weights
           # z ∈ [-1, 1] and θ ∈ [0, π] run in opposite directions, but cos(θ) == -z
           cθ = -cθ[1:nθh]
           sθΔθ = sθΔθ[1:nθh]
           rings = Vector{Ring}(undef, nθh)
           for ii in 1:nθh
               o₁ = ii           # northern-ring offset
               o₂ = nθ - ii + 1  # southern-ring offset
               offs = (o₁, o₁ == o₂ ? 0 : o₂)  # pairs, except if both aligned with equator
               rings[ii] = Ring(offs, stride, nϕ, cθ[ii], ϕ_π, sθΔθ[ii] * Δϕ)
           end
           return rings
       end
rings (generic function with 3 methods)
```

With just `rings` implemented, generic methods will then calculate other necessary
parameters by inspecting the ring list.
For example, `synthesize` allocates and fills a map buffer, and this is handled generically
by counting the number of pixels described by the ring structures and allocating a vector
of equal length:
```jldoctest examplepix
julia> glpix = GLPixelization(100, 200)
100-ring GLPixelization{Float64} with 20000 pixels

julia> npix(glpix)
20000

julia> buffer(glpix) |> println∘summary
20000-element Vector{Float64}
```

While convenient (because only the `ring` method must be written), the generic fallback
methods are less efficient and/or produce less desirable outputs than if specialized
methods are written.

For example, the Gauss-Legendre grid is very similar to the ECP grid in that it is often
very convenient to represent it as a matrix rather than a vector.
The _shape_ of the map buffers is immaterial, as long as the offsets and strides returned
by `rings` correspond to elements within the array's linear indices.
(Or said differently, the transforms effectively see any map buffer as its equivalent
`reshape(mapbuf, :)` vector representation.)
Overloading the [`buffer`](@ref) function permits synthesizing matrix maps arrays instead
of vectors:
```jldoctest examplepix
julia> function buffer(glpix::GLPixelization, ::Type{T}) where {T}
           return Matrix{T}(undef, glpix.nθ, glpix.nϕ)
       end
buffer (generic function with 4 methods)

julia> buffer(glpix) |> println∘summary
100×200 Matrix{Float64}
```

!!! note
    Be sure to overload the 2-argument method which explicitly takes the descired element
    type.
    The 1-argument method will forward the pixelization's type parameter as the second
    argument if necessary.

The other reason to implement specializations of the
[optional interface](@ref pixelinterface)
is that generic methods infer the necessary pixelization properties by iterating through the
list of all rings and extracting/calculating the relevant quantites.
We can see this quite obviously as memory allocations when running these functions,
despite the grid size alone provides all the necessary information (and therefore should
not require actually allocating the rings vector):
```jldoctest examplepix
julia> (nring(glpix), npix(glpix), nϕmax(glpix)); # warmup

julia> @time (nring(glpix), npix(glpix), nϕmax(glpix))
  0.000063 seconds (332 allocations: 34.312 KiB)
(100, 20000, 200)
```
By overloading specializations for the `GLPixelization` type,
```jldoctest examplepix
julia> nring(glpix::GLPixelization) = glpix.nθ
nring (generic function with 3 methods)

julia> npix(glpix::GLPixelization) = glpix.nθ * glpix.nϕ
npix (generic function with 3 methods)

julia> nϕmax(glpix::GLPixelization) = glpix.nϕ
nϕmax (generic function with 3 methods)
```
we can avoid such unnecessary computation and allocation:
```jldoctest examplepix
julia> (nring(glpix), npix(glpix), nϕmax(glpix)); # warmup

julia> @time (npix(glpix), nϕmax(glpix))
  0.000009 seconds (2 allocations: 48 bytes)
(20000, 200)
```
(where the remaining allocation is due to the returned tuple of values).
