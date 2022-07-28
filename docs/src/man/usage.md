# [Usage](@id usage)
```@meta
CurrentModule = SphericalHarmonicTransforms
```

```@contents
Pages = ["usage.md"]
Depth = 2
```

## Getting Started

This library is **not** registered in Julia's
[General registry](https://github.com/JuliaRegistries/General), so the package must be
installed either by cloning it directly:
```julia-repl
(@v1.7) pkg> add https://github.com/jmert/SphericalHarmonicTransforms.jl
```
or by making use of my [personal registry](https://github.com/jmert/Registry.jl)
```julia-repl
(@v1.7) pkg> registry add https://github.com/jmert/Registry.jl

(@v1.7) pkg> add SphericalHarmonicTransforms
```
After installing, just load like any other Julia package:
```julia-repl
julia> using SphericalHarmonicTransforms
```

`SphericalHarmonicTransforms` exports only a couple of functions from its public interface:
pairs of functions for spherical harmonic transform synthesis
([`synthesize`](@ref)/[`synthesize!`](@ref)) and analysis
([`analyze`](@ref)/[`analyze!`](@ref)).

To access other portions of the [public API](../lib/public.md), the bindings must be
referred to by full name or imported into the local scope.

!!! tip
    On Julia v1.6 or newer, you can abbreviate `SphericalHarmonicTransforms` as `SHT` with
    the new `import ... as ...` syntax:
    ```julia-repl
    julia> import SphericalHarmonicTransforms as SHT
    ```

## Performing Transforms

```@setup Usage
using SphericalHarmonicTransforms
const SHT = SphericalHarmonicTransforms
```

!!! note
    See the [Developer Docs](@ref pixelnew) for more details on implementing support for
    other pixelizations.

The only built-in pixelization is for an
[equidistant cylindrical projection](https://en.wikipedia.org/wiki/Equirectangular_projection)
(ECP), represented by an instance of an [`ECPPixelization`](@ref) structure.
The constructor takes at least the number of pixels in the colatitude and azimuth
directions, and optionally an element type as a type parameter which defaults to
`Float64`.
(The type parameter is used to determine the element type of newly allocated output arrays.)

```@repl Usage
pix = SHT.ECPPixelization(50, 100) # == SHT.ECPPixelization{Float64}(50, 100)
```

The object `pix` describes an array which corresponds to a column-major ordered `Float64`
matrix of size 50 × 100 that covers the sphere with pixels that are spaced uniformly in
latitude and longitude.

With the pixelization described, a set of spherical harmonic coefficients ``a_{ℓm}`` can
be synthesized to form a map using [`synthesize`](@ref).
For example, we can draw the real part of the quadrupole moment ``a_{20}``:
```@repl Usage
alm = zeros(ComplexF64, 3, 3); alm[3, 1] = 1;  # 1-indexed array for 0-indexed quantity
quad = synthesize(pix, alm)
```
and likewise we can (approximately) analyze the rendered map back to spherical harmonic
coefficients with [`analyze`](@ref):
```@repl Usage
analyze(pix, quad, #=lmax=# 5)
```

!!! warning
    The analysis transform is defined mathematically as an integral, but the operation
    demonstrated here is a finite sum over a discrete map.
    This is the reason for seeing non-zero values in the ``a_{00}`` and ``a_{40}`` terms,
    which is due in part to the very low resolution of the grid being used.
    While higher resolutions will decrease the severity of this effect, care needs to be
    used when performing spherical harmonic transforms since analysis does not inherently
    recover input synthesis coefficients (which is notably unlike the discrete Fourier
    transform).
