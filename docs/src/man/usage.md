# [Usage](@id usage)
```@meta
CurrentModule = SphericalHarmonicTransforms
DocTestSetup = quote
    using SphericalHarmonicTransforms
    const SHT = SphericalHarmonicTransforms
end
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
```jldoctest Usage
julia> using SphericalHarmonicTransforms
```

!!! tip
    On Julia v1.6 or newer, you can abbreviate `SphericalHarmonicTransforms` as `SHT` with
    the new `import ... as ...` syntax:
    ```julia-repl
    julia> import SphericalHarmonicTransforms as SHT
    ```

## Working with pixelizations

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

```jldoctest Usage
julia> pix = SHT.ECPPixelization(50, 100) # == SHT.ECPPixelization{Float64}(50, 100)
50×100 ECPPixelization{Float64}
```

The object `pix` describes an array which corresponds to a column-major ordered `Float64`
matrix of size 50 × 100 that covers the sphere with pixels that are spaced uniformly in
latitude and longitude.
