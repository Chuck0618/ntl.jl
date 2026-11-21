using NTL
using Random
using Test
using Statistics

# This tour exercises every exported interface in the current module.

# Prime-field metadata and arithmetic -----------------------------------------
@assert zzp_Min == 2
@assert zzp_Max > zzp_Min
F = zzp(101)
@assert modulus(F) == 101
@assert characteristic(F) == 101
@assert parent(F(1)) == F
@assert info(zzp).minimum_modulus == 2
@assert info(F).modulus == 101
@assert info(F(3)).representative == 3

x, y = F(3), F(4)
@assert rep(x) == 3
@assert Int(x) == 3
@assert BigInt(x) == 3
@assert zero(F) == F(0)
@assert one(F) == F(1)
@assert oneunit(F) == F(1)
@assert zero(x) == zero(F)
@assert one(x) == one(F)
@assert oneunit(x) == one(F)
@assert iszero(zero(F))
@assert isone(one(F))
@assert !iszero(x)
@assert !isone(y)
@assert x + y == F(7)
@assert x - y == F(100)
@assert -x == F(98)
@assert x * y == F(12)
@assert x / y == F(3 * inv(y))
@assert add(x, y) == x + y
@assert sub(x, y) == x - y
@assert negate(x) == -x
@assert mul(x, y) == x * y
@assert sqr(x) == x * x
@assert divide(x, y) == x / y
@assert power(x, 3) == x^3
@assert inv(x) * x == one(F)
@assert muladd(x, y, x) == x * y + x
@assert swap(x, y) == (y, x)

# Low-level modular helpers ----------------------------------------------------
@assert AddMod(100, 2, 101) == 1
@assert SubMod(1, 2, 101) == 100
@assert NegateMod(1, 101) == 100
@assert MulMod(3, 4, 101) == 12
@assert SqrMod(3, 101) == 9
@assert InvMod(3, 101) == 34
@assert InvModStatus(3, 101) == (true, 34)
@assert InvModStatus(0, 101) == (false, 0)
@assert PowerMod(3, -1, 101) == 34
@assert PowerMod(3, big(1000), 101) == 1

@assert VectorConv(F, [-1, 0, 2, 103]) == F.([100, 0, 2, 2])
@assert length(VectorRandom(F, 4; rng = MersenneTwister(1))) == 4
@assert random_zzp(F; rng = MersenneTwister(2)) isa F
@assert VectorMulMod([1, 2, 3], 3, 101) == [3, 6, 9]
out = zeros(Int, 3)
@assert VectorMulMod!(out, [1, 2, 3], 3, 101) === out
@assert out == [3, 6, 9]

# Arbitrary-precision modular rings -------------------------------------------
R = ZZP(12)
@assert R isa ZZPField
@assert R(5) isa ZZp
@assert ZZ_p === ZZp # Legacy alias retained for compatibility.
@assert ZZp(R, 5) == R(5)
@assert modulus(R) == 12
@assert characteristic(R) == 12
@assert R(5) * R(5) == R(1)
@assert rep(R(-1)) == 11
@assert BigInt(R(5)) == 5
@assert parent(R(5)) == R
@assert zero(R) == R(0)
@assert one(R) == R(1)
@assert oneunit(R) == R(1)
@assert R(5) * inv(R(5)) == one(R)
@assert R(5)^(-1) == inv(R(5))
@test_throws DomainError inv(R(6))
@test_throws ErrorException R(1) + zzp(7)(1)
@assert length(VectorRandom(R, 2; rng = MersenneTwister(3))) == 2
@assert VectorConv(R, [-1, 2, 13]) == R.([-1, 2, 13])
@assert random_zzp(R; rng = MersenneTwister(4)) isa ZZp

# Number theory ----------------------------------------------------------------
@assert ProbPrime(101)
@assert !ProbPrime(341)
@assert MillerWitness(341, 2)
@assert !MillerWitness(101, 2)
@assert NextPrime(100) == 101
@assert PrevPrime(100) == 97
@assert ndigits(GenPrime(16; rng = MersenneTwister(5)), base = 2) == 16
@assert 0 <= RandomBits(12; rng = MersenneTwister(6)) < (big(1) << 12)
@assert ndigits(RandomLen(12; rng = MersenneTwister(7)), base = 2) == 12
@assert 0 <= RandomBnd(100; rng = MersenneTwister(8)) < 100
@assert FactorInt(99991 * 10007) == [(10007, 1), (99991, 1)]
@assert EulerPhi(36) == 12
@assert CRT([2, 3, 2], [3, 5, 7]) == (23, 105)
r = SqrRootMod(10, 13)
@assert r !== nothing && r^2 % 13 == 10
@assert SqrRootMod(2, 13) === nothing
@assert DiscreteLog(2, 8, 13) == 3
@assert PrimitiveRoot(97) == 5
@assert primitive_root(F) != 0

# Dense polynomials -------------------------------------------------------------
P = F
x = xpoly(P)
f = poly(P, [1, 2, 3])
@assert poly_mul(f, x + 1) == f * (x + 1)
@assert field(f) == P
@assert coefficients(f) == P.([1, 2, 3])
@assert coeff(f, 0) == 1
@assert coeff(f, 2) == 3
@assert coeff(f, 20) == 0
@assert degree(f) == 2 == deg(f)
@assert LeadCoeff(f) == 3
@assert ConstTerm(f) == 1
@assert leading(f) == 3
setcoeff!(f, 3, 4)
@assert coeff(f, 3) == 4
f = poly(P, [1, 2, 3])
@assert monic(2x^2 + 1) == x^2 + P(51)
@assert ismonic(monic(f))
@assert LeftShift(f, 2) == x^2 * f
@assert RightShift(LeftShift(f, 2), 2) == f
q, r = DivRem(f^2, x + 1)
@assert q * (x + 1) + r == f^2
@assert GCD((x + 1)^2, x + 1) == x + 1
d, s, t = XGCD(x^2 + 1, x + 1)
@assert d == s * (x^2 + 1) + t * (x + 1)
@assert derivative(f) == 2 + 6x
@assert diff(f) == derivative(f)
@assert evaluate(f, 2) == 17
@assert evaluate(f, [0, 1, 2]) == P.([1, 6, 17])
@assert compose(f, x + 1) == 3x^2 + 8x + 6
@assert CompMod(f, x + 1, x^2 + 1) == rem(compose(f, x + 1), x^2 + 1)
@assert truncate(f^3, 2) == poly(P, [1, 6])
@assert MulTrunc(f, x + 1, 2) == truncate(f * (x + 1), 2)
@assert truncate((1 + x + x^2) * InvTrunc(1 + x + x^2, 3), 3) == one(f)
@assert MulMod(f, x + 1, x^2 + 1) == rem(f * (x + 1), x^2 + 1)
@assert rem(InvMod(x + 1, x^2 + 1) * (x + 1), x^2 + 1) == one(f)
@assert InvModStatus(x + 1, x^2 + 1)[1]
@assert PowerMod(x + 1, 3, x^2 + 1) == rem((x + 1)^3, x^2 + 1)
@assert random_poly(P, 4) isa ModPoly
@assert BuildFromRoots(P, [1, 2]) == x^2 - 3x + 2
@assert evaluate(interpolate(P, [0, 1, 2], [1, 2, 5]), [0, 1, 2]) == P.([1, 2, 5])
@assert !isempty(SquareFreeDecomp((x + 1)^2 * (x + 2)))
@assert !isempty(DistinctDegreeDecomp((x + 1) * (x + 2)))
@assert !isempty(EqualDegreeSplit((x + 1) * (x + 2), 1))
fac = factor((x + 1)^2 * (x + 2))
@assert fac isa PolyFactorization
@assert reconstruct(fac) == (x + 1)^2 * (x + 2)
@assert all(is_irreducible(g) for (g, _) in fac.factors)
@assert is_irreducible(BuildIrred(P, 3))
@assert !isempty(FindRoots((x - 1) * (x - 2)))
@assert FindRoot(x - 2) == P(2)

# GF2 and bit-packed GF2X -------------------------------------------------------
@assert GF2(0) == GF2(0)
@assert GF2(1) + GF2(1) == GF2(0)
@assert GF2(1) * GF2(1) == GF2(1)
@assert Int(GF2(1)) == 1
@assert zero(GF2) == GF2(0)
@assert one(GF2) == GF2(1)
@assert oneunit(GF2) == GF2(1)
@assert iszero(GF2(0))
@assert isone(GF2(1))
g = GF2X(2)
@assert bits(g) == 2
@assert weight(g^2 + g + 1) == 3
@assert NumBits(g^2 + 1) == 3
@assert NumBytes(g^2 + 1) == 1
@assert random_GF2X(5) isa GF2X
@assert to_modpoly(g^2 + 1) isa ModPoly
@assert from_modpoly(to_modpoly(g^2 + 1)) == g^2 + 1
@assert BuildIrredGF2X(3) isa GF2X
@assert is_irreducible(BuildIrredGF2X(3))
@assert (g + 1)^2 == g^2 + 1
facg = factor((g + 1)^3 * (g^2 + g + 1)^2)
@assert facg isa GF2XFactorization
@assert reconstruct(facg) == (g + 1)^3 * (g^2 + g + 1)^2

# Extension fields --------------------------------------------------------------
E = extension_field(P, x^2 + 2)
@assert E isa ExtensionField
@assert ExtensionElem === ZZpE # Descriptive and NTL-facing names for the same type.
@assert ZZpE(E, gen(E)) == gen(E) # Retain the NTL-facing constructor spelling.
@assert base_field(E) == P
@assert defining_polynomial(E) == x^2 + 2
@assert extension_degree(E) == 2
@assert cardinality(E) == modulus(P)^2
@assert gen(E)^2 == E(-2)
@assert trace(gen(E)) == P(0)
@assert norm(gen(E)) == P(2)
@assert parent(gen(E)) == E
@assert rep(gen(E)) isa ModPoly
@assert modulus(gen(E)) == defining_polynomial(E)
@assert characteristic(gen(E)) == 101
@assert iszero(zero(E))
@assert isone(one(E))
@assert zero(gen(E)) == zero(E)
@assert one(gen(E)) == one(E)
@assert oneunit(gen(E)) == one(E)
@test_throws DomainError inv(zero(E))
@assert ExtensionElem(E, [1, 1]) == E([1, 1])
@assert ExtensionElem(E, 2) == E(2)
@assert ExtensionElem(E, x^3 + 1) == E(x^3 + 1)
@assert degree(rep(ExtensionElem(E, x^3 + 1))) < extension_degree(E)
@assert ExtensionElem(E, gen(E)) == gen(E)
@assert ExtensionElem(E, gen(E)) isa ExtensionElem

a = E([1, 2]) # 1 + 2t, where t = gen(E) and t^2 = -2
b = E([2, 1])
@assert coefficients(rep(a)) == P.([1, 2])
@assert a + b == E([3, 3])
@assert a - b == E([-1, 1])
@assert -a == E([-1, -2])
@assert a * b == E(rem(rep(a) * rep(b), defining_polynomial(E)))
@assert sqr(a) == a^2
@assert a * inv(a) == one(E)
@assert a / b == a * inv(b)
@assert power(a, 3) == a^3
@assert a^(-2) == inv(a)^2
@assert rand(MersenneTwister(13), E) isa ExtensionElem
@assert length(rand(MersenneTwister(14), E, 2)) == 2
u = xpoly(E)
@assert evaluate(u^2 + a, zero(E)) == a

g2 = GF2X(2)
E2 = GF2E(g2^2 + g2 + 1)
@assert cardinality(E2) == 4
@assert gen(E2)^3 == one(E2)

# Zech extension representation --------------------------------------------------
Z = zech_field(E)
@assert Z isa ZechField
z = ZechElem(Z, Z._generator)
@assert parent(z) == Z
@assert rep(z) isa Int
@assert Int(z) == rep(z)
@assert modulus(z) == defining_polynomial(E)
@assert characteristic(z) == 101
@assert extension_degree(Z) == 2
@assert cardinality(Z) == cardinality(E)
@assert extension_field(Z) == E
@assert extension_element(z) isa ExtensionElem
@assert length(Z) == cardinality(Z)
@assert length(collect(Z)) == cardinality(Z)
@assert iszero(zero(Z))
@assert isone(one(Z))
@assert zero(z) == zero(Z)
@assert one(z) == one(Z)
@assert oneunit(z) == one(Z)
@assert z * inv(z) == one(Z)
@assert z^(cardinality(Z) - 1) == one(Z)
@assert muladd(z, z, one(z)) == z * z + one(z)
@assert iszero(ZechElem(Z, 0))
@assert isone(ZechElem(Z, 1))
@assert hash(z) == hash(z)
@assert rand(MersenneTwister(9), Z) isa ZechElem
@assert rand(MersenneTwister(10), Z, 2) isa Vector

# Matrices ----------------------------------------------------------------------
A = P.([2 3 1; 4 1 5; 3 7 2])
b = P.([1, 2, 3])
@assert det(A) != 0
@assert det_exact(BigInt.(A)) != 0
@assert rank(A) == 3
@assert A * solve(A, b) == b
@assert A * matrix_inverse(A) == P.([1 0 0; 0 1 0; 0 0 1])
S = P.([1 2 3; 2 4 6; 0 1 1])
@assert iszero(det(S))
@assert rank(S) == 2
R, pivots = rref(S)
@assert pivots == [1, 2]
K = nullspace(S)
@assert length(K) == 1
@assert all(iszero, S * K[1])

B = BigInt[1 1 1; -1 0 2; 3 5 6]
result = lll(B)
@assert result.transform * B == result.basis
@assert abs(det_exact(result.transform)) == 1
@assert result.verified
@assert result.mode == :fast
@assert LLL(B).basis == result.basis
@assert gram_schmidt_exact(result.basis) isa Tuple

# NTT, generalized FFT, and acceleration -----------------------------------------
NTT_F = zzp(193) # 64 divides 193 - 1, so both 32- and 64-point transforms are valid.
omega = root_of_unity(32, NTT_F)
@assert omega^32 == one(omega)
va = rand(MersenneTwister(11), NTT_F, 32)
@assert intt(ntt(va)) == va
vb = rand(MersenneTwister(12), NTT_F, 32)
@assert ntt_convolve(va, vb) == convolve(va, vb)
@assert fft_convolve(va, vb) == convolve(va, vb)
@assert !isempty(FFT_PRIMES)

set_polynomial_acceleration!(true)
@assert polynomial_acceleration()
set_polynomial_acceleration!(false)
@assert !polynomial_acceleration()
set_polynomial_acceleration!(true)
println("API tour completed successfully.")
