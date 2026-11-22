using NTL
using Random

F = PrimeField(97)
x = xpoly(F)
f = (x + 1)^3 * (x^2 + x + 1) * (x^4 + 3x + 5)

fac = factor(f; rng = MersenneTwister(2026))
println("unit = ", fac.unit)
for (factor, multiplicity) in fac.factors
    println("degree ", degree(factor), ", multiplicity ", multiplicity, ": ", factor)
    @assert is_irreducible(factor)
end
@assert reconstruct(fac) == f

q, r = DivRem(f, x^2 + 1)
@assert q * (x^2 + 1) + r == f

m = BuildIrred(F, 6; rng = MersenneTwister(7))
@assert degree(m) == 6 && is_irreducible(m)
a = PowerMod(x + 1, big(10)^100, m)
@assert a == PowerMod(x + 1, big(10)^100, m)

xs = F.([0, 1, 2, 3])
ys = F.([1, 4, 15, 40])
p = interpolate(F, xs, ys)
@assert evaluate(p, xs) == ys
println("interpolating polynomial = ", p)
println("Polynomial example passed.")
