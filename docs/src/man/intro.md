# [Introduction](@id intro)

```@contents
Pages = ["intro.md"]
Depth = 2
```

## [Definition](@id spharms-defn)

The spherical harmonics ``Y_\ell^m(\theta,\phi)`` arise from solving Laplace's equation
```math
\begin{align}
    \nabla^2 \psi = 0
\end{align}
```
in spherical coordinates.
The equation is separable into a radial component ``R(r)`` and an angular part
``Y(\theta,\phi)`` such that the total solution is
``\psi(r,\theta,\phi) \equiv R(r)Y(\theta,\phi)``.
We ignore the radial component and continue with only the angular part.

The second-order differential equation for the angular component, written in the standard
physicists' form where ``\theta`` is the colatitude and ``\phi`` the azimuth angles, is
```math
\begin{align}
    \frac{1}{\sin\theta} \frac{\partial}{\partial\theta}
        \left( \sin\theta \frac{\partial}{\partial\theta} Y \right)
        + \frac{1}{\sin^2\theta} \frac{\partial^2}{\partial\phi^2} Y
        + \ell(\ell+1) Y &= 0
\end{align}
```
The differential equation is further separable among the two angular coordinates, with the
result being a family of solutions parameterized by the integer constants ``\ell > 0``
and ``|m| < \ell``.
Including normalization factors[^norm], the spherical harmonics
``Y_{\ell m}(\theta,\phi)`` can be written as
```math
\begin{align}
    Y_{\ell m}(\theta, \phi) &\equiv
        \sqrt{\frac{2\ell+1}{4\pi} \frac{(\ell-m)!}{(\ell+m)!}}
        P_\ell^m(\cos\theta) \, e^{im\phi}
    = \lambda_\ell^m(\cos\theta) \, e^{im\phi}
\end{align}
```
where ``\lambda_\ell^m(\cos\theta)`` are pre-normalized Associated Legendre Polynomials
for spherical harmonic transforms (see
[`AssociatedLegendrePolynomials.jl`](https://github.com/jmert/AssociatedLegendrePolynomials.jl))
and ``e^{im\phi}`` are the complex exponentials.

[^norm]:
    Be aware that there are multiple conventions for how to define the spherical harmonics
    that differ in the details of the normalization factors.
    For instance, the so-called
    [Condon-Shortley phase factor](https://mathworld.wolfram.com/Condon-ShortleyPhase.html)
    of ``(-1)^m`` can either be omitted entirely, included as an explicit term in the
    normalization of the spherical harmonics, or included as an explicit symmetry in the
    definition of the Legendre polynomials.
    We use this last option here, which is also consistent with typical physicists'
    experience in e.g. *Classical Electrodynamics* by [J.D. Jackson](@ref bib-spharms).
    See also [Abramowitz & Stegun](@ref bib-spharms).

The spherical harmonics are a complete, orthonormal basis for functions on the sphere
``(\theta,\phi) \in \mathbb{S} = [0,\pi]\times[0,2\pi]``.
Therefore they satisfy the condition that
```math
\begin{align*}
    \int_\mathbb{S} Y_{\ell m}(\theta,\phi) Y_{\ell' m'}(\theta,\phi)
        \,d\Omega &= \delta_{\ell\ell'} \delta_{mm'}
\end{align*}
```
where ``d\Omega = \sin\theta \,d\theta\,d\phi``.

It then follows that an arbitrary function on the sphere ``f(\theta,\phi)`` can be
harmonically decomposed into an ensemble of coefficients ``a_{\ell m}``:
```math
\begin{align}
    a_{\ell m} &= \int_\mathbb{S} f(\theta, \phi) \overline{Y_{\ell m}}(\theta,\phi) \,d\Omega
\end{align}
```
where the overline on ``\overline{Y_{\ell m}}`` denotes complex conjugation.
The function is (re)synthesized from the ``a_{\ell m}``s as the linear combination of the
spherical harmonics:
```math
\begin{align}
    f(\theta,\phi) &= \sum_{\ell=0}^{\ell_\mathrm{max}}
        \sum_{m=-\ell}^{\ell} a_{\ell m} Y_{\ell m}(\theta,\phi)
        \label{eqn:complex_synthesis}
\end{align}
```

## Real Symmetries

In general, ``f`` may be a complex function (i.e. ``f : \mathbb{S} \rightarrow \mathbb{C}``),
but if it is a real function (``f : \mathbb{S} \rightarrow \mathbb{R}``), then there are
some additional symmetries of note.

First, due to an underlying property of the Associated Legendre polynomials, the
negative orders (``m < 0``) are related to the positive orders according to
```math
\begin{align}
    Y_{\ell (-m)}(\theta,\phi) &= (-1)^m \,\overline{Y_{\ell (+m)}}(\theta,\phi)
\end{align}
```
Given this relationship between the positive and negative orders, the harmonic coefficients
have a similar property as well, with
```math
\begin{align}
    a_{\ell(-m)} &= (-1)^m \,\overline{a_{\ell(+m)}}
\end{align}
```
This directly constraints the order ``m = 0`` modes to being real coefficients since
``x = \overline x`` if and only if ``x`` is real.
Furthermore, this permits the synthesis to be performed using a sum over the range of
integers ``m = [0, \ell]`` instead of ``[-\ell, \ell]`` because by pairwise summing
the postive and negative orders with each other, we find:
```math
\begin{align*}
    a_{\ell(-m)} Y_{\ell (-m)} + a_{\ell(+m)} Y_{\ell (+m)} &=
        (-1)^m \,\overline{a_{\ell(+m)}} \cdot (-1)^m \,\overline{Y_{\ell (+m)}} +
        a_{\ell(+m)} Y_{\ell (+m)} \\
    {} &= \overline{a_{\ell(+m)} Y_{\ell (+m)}} + a_{\ell(+m)} Y_{\ell (+m)} \\
    {} &= 2 \operatorname{Re}\left[ a_{\ell(+m)} Y_{\ell (+m)} \right]
\end{align*}
```
Note that this excludes the ``m = 0`` cases which have no pair so there needs to be an extra
multiplicity factor in the summation — for instance, choosing a notation that uses the
Kronecker delta,
``(2-\delta_{m0}) = \begin{cases} 1 & \text{if } m=0 \\ 2 & \text{otherwise} \end{cases}``,
and
```math
\begin{align}
    f(\theta,\phi) &= \sum_{\ell=0}^{\ell_\mathrm{max}}
        \sum_{m=0}^{\ell} (2 - \delta_{m0}) a_{\ell m} Y_{\ell m}(\theta,\phi)
\end{align}
```
A similar argument applies to harmonic analysis of a field into the ``a_{\ell m}``
coefficients except there is no need for an extra multiplicity factor; instead you simply
run the integration for only the non-negative orders.

The second useful symmetry is that the Legendre polynomials are even/odd for even/odd
combinations of orders and degrees, so the spherical harmonics inherit the even/odd
symmetry across the equator:
```math
\begin{align}
    P_\ell^m(-z) = (-1)^{\ell+m} P_\ell^m(z)
    \quad\Rightarrow\quad
    Y_{\ell m}(\pi - \theta, \phi) = (-1)^{\ell+m} Y_{\ell m}(\theta,\phi)
\end{align}
```
This is a useful property when it comes to designing pixelization schemes on the sphere
since symmetry over the equator allows for reuse in the calculation of Legendre polynomials
for pairs of latitudes with only a much [computationally] cheaper negation required of
one of half of the pair.

## FFT-based Ring Transforms

If the sphere is discretized by a series of isolatitude rings, each of which are comprised
of uniformly-spaced pixels in azimuth, then the fast Fourier transform (FFT) may be
employed to reduce the computational complexity of the transform.

Starting from Eqn. ``\ref{eqn:complex_synthesis}``, changing the order of summation so
that the sums over degrees ``\ell`` are done before the order ``m`` reveals the inner
summation in a form which corresponds to a discrete Fourier transform.
```math
\begin{align*}
    \def\ellmax{\ell_{\mathrm{max}}}
    f(\theta,\phi) &= \sum_{m=-\ellmax}^{\ellmax}
        \left[ \sum_{\ell=|m|}^{\ellmax} a_{\ell m} \lambda_\ell^m(\cos\theta) \right]
        e^{im\phi}
\end{align*}
```
Then for uniformly spaced pixels in the ring, the azimuth coordinates are
```math
\begin{align*}
    \phi_k = \phi_0 + \frac{2\pi k}{N_\phi}
\end{align*}
```
where ``\phi_0`` is the offset of the first pixel in the ring away from ``\phi = 0``,
``N_\phi`` is the number of pixels in the ring, and ``k`` is an integer in
``[0, N_\phi - 1]``.
Combined with the reordered summation, the expression can be rewritten as
```math
\begin{align*}
    \def\ellmax{\ell_{\mathrm{max}}}
    f(\theta,\phi_k) &= \sum_{m=-\ellmax}^{\ellmax} g_m e^{2\pi imk/N_\phi}
    &\text{where }
        g_m \equiv e^{im\phi_0} \sum_{\ell=|m|}^{\ellmax} a_{\ell m} \lambda_\ell^m(\cos\theta)
\end{align*}
```
The summation on the left is discrete Fourier transform of the coefficients ``g_m``.
When ``N_\phi = 2\ellmax + 1``, the correspondance to the FFT is exact, and
``f(\theta, \phi_k) = \mathcal{F}^{-1}(g_m)``, using
[FFTW](http://www.fftw.org/)'s sign convention and ``\mathcal{F}`` and ``\mathcal{F}^{-1}``
are the foward and inverse FFT operations, respectively.
Also recall that "negative frequencies" are just a matter of interpretation and are
equivalent to a higher positive frequency one period later.
For instance, the ``m = -1`` term is equivalent in the Fourier transform to the
``(N_\phi - 1)``-th Fourier frequency.

For cases when ``N_\phi \neq 2\ellmax + 1``, the summation must be properly zero-padded
or alias-truncated — for further details, see
[Notes on Calculating the Spherical Harmonics (Spherical Harmonics Series, Part I)](@ref bib-spharms).
