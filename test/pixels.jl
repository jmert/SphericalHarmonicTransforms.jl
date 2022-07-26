import ..TestSuite
const SHT = SphericalHarmonicTransforms

# even & odd sizes
TestSuite.runtests(SHT.ECPPixelization(10, 20))
TestSuite.runtests(SHT.ECPPixelization(11, 21))
