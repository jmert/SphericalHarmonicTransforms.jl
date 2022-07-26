using SphericalHarmonicTransforms
import ..TestSuite: synthesize_reference, analyze_reference

# Generate complex fields from running real-only synthesis on real and imaginary
# components separately.
function synthesize_reference_complex(ℓ, m, θ, ϕ, T::Type = Float64)
    alms = zeros(Complex{T}, ℓ + 1, m + 1)
    # Assign unit power to single delta (ℓ,m) mode.
    #
    # Real part
    alms[ℓ+1,m+1] = 1
    ref = synthesize_reference(alms, θ, ϕ)
    # Complex part -- Use (-im) not (+im) since the goal is to swap the internal
    #                 (a + ib) to (b + ia), which requires multiplication by -im.
    alms[ℓ+1,m+1] *= -im
    ref += synthesize_reference(alms, θ, ϕ) .* im
    # account for doubling due to assumption of real-only symmetry in synthesize_reference
    ref ./= m == 0 ? 1 : 2
    return ref
end

# Because of the disparity between analysis as integration of a field and summation over
# a discrete map, a single round of analyzing any given map suffers from non-floating
# point errors. By iterating and repeatedly analyzing the residual between input and
# analyze+(re)synthesize, the harmonic coefficient slowly converge to what we'd get if
# we had much higher resolution input maps.
#
# Generically, this is a problem that falls under the umbrella of "quadrature weights", but
# that's a complex topic that we bypass by just brute-force converging our solution.
function analyze_reference_complex_iter(map, θ, ϕ, lmax, mmax)
    R = eltype(θ)
    tol = sqrt(eps(maximum(abs, map)))
    # norm is the correct normalization for an ECP pixelization
    norm = 2R(π)^2 / length(θ)

    function iter(ecp)
        local alms
        alms = norm .* analyze_reference(ecp, θ, ϕ, lmax, mmax)
        for ii in 2:10  # allow up to 10 iterations before bailing
            ecp′ = synthesize_reference(alms, θ, ϕ)
            δecp = ecp - ecp′
            # check for convergence
            if maximum(abs, δecp) < tol
                break
            end
            alms .+= norm .* analyze_reference(δecp, θ, ϕ, lmax, mmax)
        end
        return alms
    end
    # analyze each of real and imaginary parts separately since the routines require
    # real inputs.
    alms = iter(real.(map)) .+ im .* iter(imag.(map))
    return alms
end

# Define the analytic expressions for the first few spherical harmonics; verify
# implementations against these.
Y00(θ, ϕ) =  sqrt(1 / 4oftype(θ, π)) * complex(true)
Y10(θ, ϕ) =  sqrt(3 / 4oftype(θ, π)) * cos(θ) * complex(true)
Y11(θ, ϕ) = -sqrt(3 / 8oftype(θ, π)) * sin(θ) * cis(ϕ)
Y20(θ, ϕ) =  sqrt(5 / 16oftype(θ, π)) * (3cos(θ)^2 - 1) * complex(true)
Y21(θ, ϕ) = -sqrt(15 / 8oftype(θ, π)) * sin(θ) * cos(θ) * cis(ϕ)
Y22(θ, ϕ) =  sqrt(15 / 32oftype(θ, π)) * sin(θ)^2 * cis(2ϕ)

# Use simple ECP grid for tests
ecprange(l, h, n) = (Δ = (h - l) / n; range(Δ / 2, h - Δ / 2, length = n))
Npix = 50
θ = repeat(ecprange(0.0, 1.0π, Npix), 1, 2Npix)
ϕ = repeat(ecprange(0.0, 2.0π, 2Npix)', Npix, 1)

# semi-arbitrary choice of limits, chosen only to be higher than Y22 as a mild check
# against obvious problems
lmax, mmax = 10, 3


@testset "Analytic checks — synthesis" begin
    # Validates the baseline reference implementation
    @test Y00.(θ, ϕ) ≈ synthesize_reference_complex(0, 0, θ, ϕ)
    @test Y10.(θ, ϕ) ≈ synthesize_reference_complex(1, 0, θ, ϕ)
    @test Y11.(θ, ϕ) ≈ synthesize_reference_complex(1, 1, θ, ϕ)
    @test Y20.(θ, ϕ) ≈ synthesize_reference_complex(2, 0, θ, ϕ)
    @test Y21.(θ, ϕ) ≈ synthesize_reference_complex(2, 1, θ, ϕ)
    @test Y22.(θ, ϕ) ≈ synthesize_reference_complex(2, 2, θ, ϕ)
end

@testset "Analytic checks — analysis" begin
    # Validates the baseline reference implementation
    analyze(Y) = analyze_reference_complex_iter(Y.(θ, ϕ), θ, ϕ, lmax, mmax)
    expect(ℓ, m) = setindex!(zeros(complex(eltype(θ)), lmax+1, mmax+1), 1.0, ℓ+1, m+1)
    @test analyze(Y00) ≈ expect(0, 0)
    @test analyze(Y10) ≈ expect(1, 0)
    @test analyze(Y11) ≈ expect(1, 1)
    @test analyze(Y20) ≈ expect(2, 0)
    @test analyze(Y21) ≈ expect(2, 1)
    @test analyze(Y22) ≈ expect(2, 2)
end
