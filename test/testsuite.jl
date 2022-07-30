module TestSuite
    using Test, SphericalHarmonicTransforms
    import SphericalHarmonicTransforms:
        AbstractRingPixelization, Ring, rings, buffer, npix, nϕmax

    include("helpers/reference.jl")

    # returns
    # - θ, ϕ = coordinates of the pixels in the pixelization
    # - Δ = the integration-factor that normalizes a map being passed to analyze_reference;
    #       this is ΔΩ / sin(θ) since analyze_reference includes the multiplication by
    #       sin(θ) as part of its implementation
    function ringpix2coords(ringpix::AbstractRingPixelization)
        θ = buffer(ringpix)
        ϕ = buffer(ringpix)
        Δ = buffer(ringpix)
        R = eltype(θ)
        for rpix in rings(ringpix)
            o₁, o₂ = rpix.offset
            N = rpix.nϕ
            cθ = rpix.cosθ
            sθ = sqrt(fma(cθ, -cθ, one(cθ)))
            _strided(o) = range(o, step = rpix.stride, length = N)
            if o₁ != 0
                @view(θ[_strided(o₁)]) .= acos(rpix.cosθ)
                @view(ϕ[_strided(o₁)]) .= R(π) .* (rpix.ϕ₀_π .+ 2 .* (0:N-1) ./ N)
                @view(Δ[_strided(o₁)]) .= rpix.ΔΩ / sθ
            end
            if o₂ != 0
                @view(θ[_strided(o₂)]) .= acos(-rpix.cosθ)
                @view(ϕ[_strided(o₂)]) .= R(π) .* (rpix.ϕ₀_π .+ 2 .* (0:N-1) ./ N)
                @view(Δ[_strided(o₂)]) .= rpix.ΔΩ / sθ
            end
        end
        return (θ, ϕ, Δ)
    end

    function runtests(ringpix; rtol = nothing)
        @testset "Test Suite - $(typeof(ringpix))" begin
            interface(ringpix)
            reference(ringpix, rtol)
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

    function reference(ringpix::AbstractRingPixelization{R}, rtol) where {R}
        if rtol === nothing
            rtol = sqrt(eps(one(R)))
        end
        @testset "Reference Comparisons" begin
            θ, ϕ, Δ = ringpix2coords(ringpix)
            lmax = 5
            alms = zeros(complex(R), lmax + 1, lmax + 1)
            @testset "a_{$ℓ,$m}" for ℓ in 0:lmax, m in 0:ℓ
                fill!(alms, zero(R))
                alms[ℓ+1, m+1] = one(R)

                mapbuf = synthesize_reference(alms, θ, ϕ)
                alms′ = analyze_reference(mapbuf .* Δ, θ, ϕ, lmax, lmax)
                @test synthesize(ringpix, alms) ≈ mapbuf rtol=rtol
                @test analyze(ringpix, mapbuf, lmax, lmax) ≈ alms′ rtol=rtol
            end
        end
        return nothing
    end
end
