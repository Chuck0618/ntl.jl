using NTL
using Random

F = zzp(101)
rng = MersenneTwister(2026)

# Standard Julia random arrays work directly with the field element type.
v = rand(rng, F, 8)
w = VectorRandom(F, 8; rng = rng)
println("random v = ", v)
println("random w = ", w)

# VectorConv accepts ordinary and arbitrary-precision integers.
coefficients = VectorConv(F, [-1, 0, 102, big(10)^40])
println("converted coefficients = ", coefficients)

# Componentwise arithmetic uses Julia broadcasting.
@assert v + w == v .+ w
@assert 2 .* v == v .+ v
@assert all(x -> 0 <= rep(x) < modulus(F), coefficients)

# The low-level vector routine works with reduced Int representatives.
x = [0, 1, 50, 100]
y = VectorMulMod(x, 7, 101)
@assert y == [MulMod(a, 7, 101) for a in x]

VectorMulMod!(x, x, 7, 101)  # In-place aliasing is supported.
@assert x == y

println("All vector and random-generation checks passed.")
