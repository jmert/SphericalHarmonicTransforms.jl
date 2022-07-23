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
