using NTL
using Random

n = big(2)^89 - 1
@assert ProbPrime(n)
println("2^89-1 is prime: ", n)

composite = big(99991) * big(10007)
factors = FactorInt(composite)
println("factorization = ", factors)
@assert prod(p^e for (p, e) in factors; init = BigInt(1)) == composite

x, M = CRT([2, 3, 2], [3, 5, 7])
@assert (x, M) == (23, 105)
println("CRT solution = ", x, " modulo ", M)

root = SqrRootMod(10, 13)
@assert root !== nothing && root^2 % 13 == 10
println("sqrt(10) mod 13 = ", root)

g = PrimitiveRoot(97)
e = 42
y = powermod(g, e, 97)
recovered = DiscreteLog(g, y, 97)
@assert recovered !== nothing && powermod(g, recovered, 97) == y
println("discrete log = ", recovered)

p = GenPrime(48; rng = MersenneTwister(2026))
@assert ProbPrime(p) && ndigits(p, base = 2) == 48
println("generated 48-bit prime = ", p)
