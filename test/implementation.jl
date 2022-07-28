import Random
import SphericalHarmonicTransforms: ECPPixelization, buffer
import ..TestSuite: synthesize_reference, analyze_reference, ringpix2coords

"""
    gen_alms(::Type{T}, lmax::Integer, mmax::Integer = lmax; seed = nothing) where T

Generates random alms up to ``(ℓ,m) = (lmax, mmax)`` of type `complex(T)`. Note that the
DC (``Y_{00}``) term is not set.

A seed value may be passed as a keyword argument to reset the random number generator before
generating alms. Values are drawn such that the return values for ``lmax < lmax′`` given
the same seed will be equal for all ``(ℓ, m)`` where ``ℓ ≤ lmax``.
"""
function gen_alms(::Type{T}, lmax::Integer, mmax::Integer = lmax; seed = nothing) where T
    if seed !== nothing
        Random.seed!(seed)
    end
    alms = zeros(complex(T), lmax + 1, mmax + 1)
    @inbounds for l in 1:lmax
        alms[l+1,1] = randn(real(T))
        alms[l+1,2:l+1] .= randn.(complex(T))
    end
    return alms
end

# even grid sizes
Npix = 50
ecppix = ECPPixelization(Npix, 2Npix)
θ, ϕ, _ = ringpix2coords(ecppix)
# odd grid sizes
ecppix′ = ECPPixelization(Npix+1, 2Npix+1)
θ′, ϕ′, _ = ringpix2coords(ecppix′)

# two sets of alms: one that is bandlimited below map Nyquist, and one that extends well
# above the Nyquist
lmax_lo = Npix ÷ 2 - 1
lmax_hi = 2Npix
alms_lo = gen_alms(Float64, lmax_lo, seed = 5)
alms_hi = gen_alms(Float64, lmax_hi, seed = 5)

@testset "Aliased ring synthesis" begin
    # Checks that the fast FFT algorithm correctly includes the effects of aliasing
    # high-m modes back down to low-m modes
    ref = synthesize_reference(alms_hi, θ, ϕ)
    ecp = synthesize(ecppix, alms_hi)
    @test ref ≈ ecp
    # repeat again with one longer index to make sure even and odd cases both handled
    # correctly.
    ref = synthesize_reference(alms_hi, θ′, ϕ′)
    ecp = synthesize(ecppix′, alms_hi)
    @test ref ≈ ecp

    # nϕ == 1 is a special case in aliasing that must be handled appropriately
    alm00 = ComplexF64[1 0; 0 0]
    ref = synthesize_reference(alm00, θ[:,1:1], fill(π/1, Npix, 1))
    ecp = synthesize(ECPPixelization(Npix, 1), alm00)
    @test ref == ecp == fill(1/sqrt(4π), Npix, 1)
    # make sure nϕ == 2 still works
    alm11 = ComplexF64[0 0; 0 1im]
    ref = synthesize_reference(alm11, θ[:,1:2], repeat([π/2 3π/2], Npix, 1))
    ecp = synthesize(ECPPixelization(Npix, 2), alm11)
    @test ref ≈ ecp
end

@testset "Aliased ring analysis" begin
    a22 = ComplexF64[0 0 0; 0 0 0; 0 0 1]
    # Checks that the fast FFT algorithm correctly pads the rings if necessary to reach
    # analyzing to high-m.
    srcmap = synthesize(ecppix, a22)
    ref = analyze_reference(srcmap, θ, ϕ, lmax_hi) .* (2π^2 / prod(size(θ)))
    ecp = analyze(ecppix, srcmap, lmax_hi)
    @test ref ≈ ecp
    # repeat again with one longer index to make sure even and odd cases both handled
    # correctly.
    srcmap′ = synthesize(ecppix′, a22)
    ref = analyze_reference(srcmap′, θ′, ϕ′, lmax_hi) .* (2π^2 / prod(size(θ′)))
    ecp = analyze(ecppix′, srcmap′, lmax_hi)
    @test ref ≈ ecp
end

@testset "Type promotion" begin
    pix32 = ECPPixelization{Float32}(5, 10)
    alms64 = ComplexF64[1 0; 0 0]
    map64 = buffer(pix32, Float64)

    @test @inferred(synthesize(pix32, alms64)) isa Matrix{Float64}
    @test @inferred(analyze(pix32, map64, 1)) isa Matrix{ComplexF64}
    # TODO: Figure out how to test the internal type promotion
    #pix64 = ECPixelization{Float64}(5, 10)
    #@test @inspect_inner(synthesize!(map32, pix64, alms32)) == Float64
    #@test @inspect_inner(analyze!(alms32, pix64, map32, 1)) == Float64
end
