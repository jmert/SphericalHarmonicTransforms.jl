using Documenter

# Importing the following into Main changes how types are printed (and therefore impacts
# doctests) — namely, the module is not prefixed on the type which keeps the names short.
using SphericalHarmonicTransforms
using SphericalHarmonicTransforms: ECPPixelization, RingPixelization

doctest = "--fix"  in ARGS ? :fix :
          "--test" in ARGS ? true : false

DocMeta.setdocmeta!(SphericalHarmonicTransforms, :DocTestSetup, :(using SphericalHarmonicTransforms); recursive=true)

makedocs(
    format = Documenter.HTML(mathengine = Documenter.MathJax3()),
    sitename = "Spherical Harmonic Transforms",
    authors = "Justin Willmert",
    modules = [SphericalHarmonicTransforms],
    doctest = doctest,
    checkdocs = :exports,
    doctestfilters = Regex[
        r"Ptr{0x[0-9a-f]+}",
        r"[0-9\.]+ seconds \(.*\)"
    ],
    pages = [
        "SphericalHarmonicTransforms.jl Documentation" => "index.md",
        "Spherical Harmonic Transforms" => [
            "Introduction" => "man/intro.md",
            "Usage" => "man/usage.md",
            "Developer Documentation" => "man/devdocs.md",
            "Literature/References" => "man/references.md"
        ],
        "API Reference" => "lib/public.md",
    ],
)

deploydocs(
    repo = "github.com/jmert/SphericalHarmonicTransforms.jl.git",
    push_preview = true
)
