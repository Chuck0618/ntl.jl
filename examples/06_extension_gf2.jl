using NTL
using Random

# GF(2^8) using the AES polynomial x^8+x^4+x^3+x+1.
x = GF2X(2)
aes_modulus = x^8 + x^4 + x^3 + x + 1
@assert is_irreducible(aes_modulus)
E = GF2E(aes_modulus)
a = gen(E) + 1
@assert cardinality(E) == 256
@assert a^255 == one(E)
@assert a * inv(a) == one(E)
println("GF(2^8) element = ", a)
println("trace = ", trace(a), ", norm = ", norm(a))

# A prime-field extension F_5^2 = F_5[t]/(t^2+2).
F = zzp(5)
t_poly = xpoly(F)
E5 = extension_field(F, t_poly^2 + 2)
t = gen(E5)
@assert t^2 == E5(3)
@assert norm(t) == F(2)
@assert trace(t) == F(0)
println("norm(t) = ", norm(t))

# Polynomials over the extension field use the same generic algorithms.
y = xpoly(E5)
f = (y + 1)^2 * (y + t)
fac = factor(f; rng = MersenneTwister(11))
@assert reconstruct(fac) == f
println("Extension-field polynomial factorization passed.")
