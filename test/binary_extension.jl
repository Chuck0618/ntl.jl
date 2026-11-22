@testset "GF2 and bit-packed GF2X" begin
    for a in 0:1, b in 0:1
        @test rep(GF2(a) + GF2(b)) == xor(a, b)
        @test rep(GF2(a) * GF2(b)) == a & b
        @test GF2(a) - GF2(b) == GF2(a) + GF2(b)
    end
    @test inv(GF2(1)) == GF2(1)
    @test_throws DomainError inv(GF2(0))
    @test GF2(1)^100 == GF2(1)
    @test rand(MersenneTwister(1), GF2, 5) == rand(MersenneTwister(1), GF2, 5)

    x = GF2X(2)
    @test zero(GF2X) == GF2X()
    @test one(GF2X) == GF2X(1)
    @test degree(x) == 1
    @test coeff(x, 1) == GF2(1)
    @test coefficients(x + 1) == GF2.([1, 1])
    @test GF2X(GF2.([1, 1, 0, 1])) == x^3 + x + 1
    @test (x + 1)^2 == x^2 + 1
    @test weight(x^5 + x^2 + 1) == 3
    @test NumBits(x^5 + 1) == 6
    @test NumBytes(x^5 + 1) == 1
    @test LeftShift(x + 1, 3) == x^4 + x^3
    @test RightShift(x^4 + x^3 + 1, 2) == x^2 + x
    @test reverse(x^3 + x + 1) == x^3 + x^2 + 1
    @test derivative(x^4 + x^3 + x + 1) == x^2 + 1
    @test evaluate(x^3 + x + 1, 0) == GF2(1)
    @test evaluate(x^3 + x + 1, 1) == GF2(1)
    @test compose(x^2 + 1, x + 1) == (x + 1)^2 + 1
    @test truncate(x^5 + x^2 + 1, 3) == x^2 + 1
    @test MulTrunc(x^4 + x, x^2 + 1, 3) == truncate((x^4 + x) * (x^2 + 1), 3)

    f = x^7 + x^3 + 1
    g = x^3 + x + 1
    q, r = DivRem(f, g)
    @test q * g + r == f
    @test div(f, g) == q
    @test rem(f, g) == r
    @test f % g == r
    d, s, t = XGCD(g, x + 1)
    @test d == s * g + t * (x + 1)
    invf = InvMod(x + 1, g)
    @test MulMod(invf, x + 1, g) == one(f)
    @test InvModStatus(x + 1, g) == (true, invf)
    @test PowerMod(x + 1, 20, g) == rem((x + 1)^20, g)
    @test PowerMod(x + 1, -1, g) == invf
    h = x^4 + x + 1
    @test truncate(h * InvTrunc(h, 8), 8) == one(h)
    @test degree(random_GF2X(12; rng = MersenneTwister(4))) <= 11

    f = (x + 1)^3 * (x^2 + x + 1)^2
    fac = factor(f; rng = MersenneTwister(8))
    @test reconstruct(fac) == f
    @test all(is_irreducible(g) for (g, _) in fac.factors)
    @test from_modpoly(to_modpoly(f)) == f
    irreducible = BuildIrredGF2X(6; rng = MersenneTwister(9))
    @test degree(irreducible) == 6
    @test is_irreducible(irreducible)
    @test count(r -> r == GF2(1), FindRoots((x + 1)^2; rng = MersenneTwister(10))) == 2
end

@testset "Prime and binary extension fields" begin
    for F in (PrimeField(5), BigModRing(7), GF2)
        x = xpoly(F)
        m = F == PrimeField(5) ? x^2 + 2 : (F == BigModRing(7) ? x^2 + 1 : x^2 + x + 1)
        E = extension_field(F, m)
        a = gen(E)
        @test cardinality(E) == cardinality(F)^2
        @test extension_degree(E) == 2
        @test base_field(E) == F
        @test defining_polynomial(E) == m
        direct = ExtensionElem(E, m + x) # Public polynomial construction must reduce first.
        @test direct == a
        @test degree(rep(direct)) < extension_degree(E)
        @test a * inv(a) == one(E)
        @test a^0 == one(E)
        @test a^-1 == inv(a)
        @test a + 2 == E(2) + a
        @test 2a == a + a
        @test trace(a) isa typeof(zero(F))
        @test norm(a) isa typeof(zero(F))
        @test norm(a) == F(ConstTerm(m)) # monic quadratic: norm(x)=constant term
        @test rand(MersenneTwister(1), E, 3) == rand(MersenneTwister(1), E, 3)
        @test parent(rand(MersenneTwister(2), E)) == E
        @test_throws DomainError inv(zero(E))

        y = xpoly(E)
        f = (y + 1)^2 * (y + a)
        fac = factor(f; rng = MersenneTwister(12))
        @test reconstruct(fac) == f
    end

    x = GF2X(2)
    E = GF2E(x^2 + x + 1)
    a = gen(E)
    @test a^3 == one(E)
    @test trace(a) == GF2(1)
    @test norm(a) == GF2(1)
    @test_throws ArgumentError GF2E(x^2 + 1)
end
