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

```@setup examplepix
using SphericalHarmonicTransforms
using SphericalHarmonicTransforms: ECPPixelization
```

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

```@repl examplepix
import SphericalHarmonicTransforms: AbstractRingPixelization, Ring,
        rings, buffer, nring, npix, nϕmax
struct GLPixelization{R} <: AbstractRingPixelization{R}
    nθ::Int
    nϕ::Int
end
GLPixelization(nθ, nϕ) = GLPixelization{Float64}(nθ, nϕ)
```

Then, the only other **required** step is to implement a specialization of the
[`rings`](@ref) method to return a vector of [`Ring`](@ref) structures that describe the
[properties of the ring](@ref pixelproperties).

```@repl examplepix
using FastGaussQuadrature: gausslegendre
function rings(glpix::GLPixelization{R}) where {R}
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
```

With just `rings` implemented, generic methods will then calculate other necessary
parameters by inspecting the ring list.

For example, we can observe that the Gauss-Legendre grid has better properties with respect
to the finite integration when compared to the ECP grid.
```@repl examplepix
alms = zeros(ComplexF64, 6, 6); alms[3, 1] = 1;  # a_{20} = 1
ecppix, glpix = ECPPixelization(50, 100), GLPixelization(50, 100);
ecpmap, glmap = synthesize(ecppix, alms), synthesize(glpix, alms);
analyze(ecppix, ecpmap, #=lmax=# 5)
analyze(glpix, glmap, #=lmax=# 5)
```
Notice that the ECP grid analysis "leaks" power into neighboring spherical harmonic
coefficients, but the GL grid almost exactly recovers the input coefficients.

Returning to the optional components of abstract pixelization interface, compare the shapes
of the two maps created above:
```@repl examplepix
size(ecpmap)
size(glmap)
```
Like the ECP grid, it may be preferable to conceptualize and represent the GL grid as
a matrix rather than a vector.
(Fundamentally, the _shape_ of the array is immaterial to `synthesize` and `analyze`, as
long as the described ring pixelization is within bounds of the array's linear indices.)
Overloading the [`buffer`](@ref) function provides a mechanism to allocate the synthesize
map with a preferred shape:
```@repl examplepix
function buffer(glpix::GLPixelization, ::Type{T}) where {T}
    return Matrix{T}(undef, glpix.nθ, glpix.nϕ)
end
size(synthesize(glpix, alms))
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
```@repl examplepix
(nring(glpix), npix(glpix), nϕmax(glpix)); # warmup
@time (nring(glpix), npix(glpix), nϕmax(glpix))
```
By overloading specializations for the `GLPixelization` type,
```@repl examplepix
nring(glpix::GLPixelization) = glpix.nθ
npix(glpix::GLPixelization) = glpix.nθ * glpix.nϕ
nϕmax(glpix::GLPixelization) = glpix.nϕ
```
we can avoid such unnecessary computation and allocation:
```@repl examplepix
(nring(glpix), npix(glpix), nϕmax(glpix)); # warmup
@time (npix(glpix), nϕmax(glpix))
```
(where the remaining allocation is due to the returned tuple of values).
