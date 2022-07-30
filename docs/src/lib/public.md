```@meta
CurrentModule = SphericalHarmonicTransforms
```
# Public API

## Functions

```@docs
synthesize
synthesize!
```

```@docs
analyze
analyze!
```

## Ring-based Pixelizations

### Pixelizations

```@docs
RingPixelization
ECPPixelization
```

### Abstract Pixelization Interface

The following functions and types are considered part of the public interface, though
they are unlikely to be used unless you are
[implementing a new pixelization](@ref pixelnew).

```@docs
AbstractRingPixelization
Ring
rings
nring
npix
nϕmax
buffer
```
