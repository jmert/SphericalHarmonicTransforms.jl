module TestSuite
    using Test, SphericalHarmonicTransforms
    import SphericalHarmonicTransforms:
        AbstractRingPixelization, Ring, rings, buffer, npix, nϕmax

    include("helpers/reference.jl")

    function runtests(ringpix)
        @testset "Test Suite - $(typeof(ringpix))" begin
            interface(ringpix)
        end
        return nothing
    end

    function interface(ringpix)
        @testset "AbstractRingPixelization Interface" begin
            @test ringpix isa AbstractRingPixelization

            # Verify basic properties which may be overloaded against the generic implementation
            @test nϕmax(ringpix) == invoke(nϕmax, Tuple{AbstractRingPixelization}, ringpix)
            @test npix(ringpix) == invoke(npix, Tuple{AbstractRingPixelization}, ringpix)

            # buffer() returns a buffer with an overridable element type
            @test eltype(buffer(ringpix, Float64)) === Float64
            @test eltype(buffer(ringpix, Float32)) === Float32

            rr = rings(ringpix)
            # if rr is a lazy computed iterable, the result must be invariant to repeated use
            @test collect(rr) == collect(rr)
        end
        return nothing
    end
end
