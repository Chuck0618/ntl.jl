using NTL
using Random

# 257-1=256 supports both the 32-point round trip and 256-point convolution.
F = zzp(257)
omega = root_of_unity(32, F)
@assert omega^32 == one(F) && omega^16 != one(F)

rng = MersenneTwister(2026)
a = rand(rng, F, 32)
@assert intt(ntt(a)) == a
println("32-point NTT round trip passed")

u = rand(rng, F, 100)
v = rand(rng, F, 80)
w = ntt_convolve(u, v)
@assert w == convolve(u, v)
println("NTT convolution length = ", length(w))
println("First five coefficients = ", w[1:5])
