# Developer Documentation

```@meta
CurrentModule = SphericalHarmonicTransforms
```

```@contents
Pages = ["devdocs.md"]
Depth = 2
```

## [Adding New Pixelizations](@id pixelnew)

### Necessary Properties of Ring-based Pixelizations

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
