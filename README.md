# SphericalHarmonicTransforms.jl

| **Documentation**                                                         | **Build Status**                                     |
|:-------------------------------------------------------------------------:|:----------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][ci-img]][ci-url][![][codecov-img]][codecov-url] |

### Installation and usage

This library is **not** registered in Julia's [General registry][General.jl],
so the package must be installed either by cloning it directly:
```
pkg> add https://github.com/jmert/SphericalHarmonicTransforms.jl
```
or by making use of my [personal registry][Registry.jl]:
```
pkg> registry add https://github.com/jmert/Registry.jl
pkg> add SphericalHarmonicTransforms
```
After installing, just load like any other Julia package:
```julia
julia> using SphericalHarmonicTransforms
```

[General.jl]: https://github.com/JuliaRegistries/General
[Registry.jl]: https://github.com/jmert/Registry.jl

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://jmert.github.io/SphericalHarmonicTransforms.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://jmert.github.io/SphericalHarmonicTransforms.jl/dev

[ci-img]: https://github.com/jmert/SphericalHarmonicTransforms.jl/actions
[ci-url]: https://github.com/jmert/SphericalHarmonicTransforms.jl/workflows/CI/badge.svg

[codecov-img]: https://codecov.io/gh/jmert/SphericalHarmonicTransforms.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jmert/SphericalHarmonicTransforms.jl
