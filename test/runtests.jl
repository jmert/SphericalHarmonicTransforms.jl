using Test, TestSetExtensions

# Importing the following into Main changes how types are printed (and therefore impacts
# doctests) — namely, the module is not prefixed on the type which keeps the names short.
using SphericalHarmonicTransforms
using SphericalHarmonicTransforms: ECPPixelization, RingPixelization

include("testsuite.jl")

function prettytime(t)
    v, u = t < 1e3 ? (t, "ns") :
           t < 1e6 ? (t/1e3, "μs") :
           t < 1e9 ? (t/1e6, "ms") :
                     (t/1e9, "s")
    return string(round(v, digits=3), " ", u)
end
macro include(file, desc)
    mod = gensym(first(splitext(file)))
    quote
        print($desc, ": ")
        t0 = time_ns()
        @testset $desc begin
            @eval module $mod
                using Test, SphericalHarmonicTransforms
                include($file)
            end
        end
        printstyled("  ", prettytime(time_ns() - t0), "\n", color=:light_black)
    end
end

@testset ExtendedTestSet "SphericalHarmonicTransforms" begin
    @include "reference.jl" "Reference Implementation"
    @include "pixels.jl" "Pixelizations"
    @include "implementation.jl" "Implementation details"
    @include "doctests.jl" "Doctests"
end
