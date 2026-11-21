@testset "Dense polynomial arithmetic" begin
    for p in (2, 3, 5, 7)
        F = zzp(p)
        x = xpoly(F)
        @test degree(poly(F, Int[])) == -1
        @test degree(x) == 1
        @test coeff(x, 1) == 1
        @test coeff(x, 0) == 0
        @test coeff(x, 20) == 0
        @test LeadCoeff(x) == 1
        @test ConstTerm(x + 2) == 2
        @test (x + 1) * (x - 1) == x^2 - 1
        @test LeftShift(x + 1, 3) == x^4 + x^3
        @test RightShift(x^4 + x^3 + x + 1, 2) == x^2 + x
        @test derivative(x^3 + 2x + 1) == 3 * x^2 + 2
        @test diff(x^2) == derivative(x^2)
        @test evaluate(x^2 + 2x + 1, 2) == (F(2) + 1)^2
        @test evaluate(x^2 + 1, [0, 1, 2]) == F.([1, 2, 5])
        @test compose(x^2 + 1, x + 1) == (x + 1)^2 + 1
        @test CompMod(x^3 + 1, x + 1, x^2 + 1) == rem(compose(x^3 + 1, x + 1), x^2 + 1)

        f = x^5 + 3x^2 + 1
        g = x^2 + 1
        q, r = DivRem(f, g)
        @test q * g + r == f
        @test degree(r) < degree(g)
        @test div(f, g) == q
        @test rem(f, g) == r
        @test f % g == r
        @test monic(-x^2 + 1) == x^2 - 1
        @test ismonic(monic(f))
        d = gcd((x + 1)^2 * (x + 2), (x + 1) * (x + 3))
        @test d == (p == 2 ? (x + 1)^2 : x + 1)
        @test GCD(f, g) == gcd(f, g)
        modulus_poly = x^2 + x + 1
        d, s, t = XGCD(modulus_poly, x + 1)
        @test d == s * modulus_poly + t * (x + 1)
        invg = InvMod(x + 1, modulus_poly)
        @test MulMod(invg, x + 1, modulus_poly) == one(f)
        ok, value = InvModStatus(x + 1, modulus_poly)
        @test ok && value == invg
        @test !InvModStatus(modulus_poly, modulus_poly)[1]
        @test PowerMod(x + 1, 10, modulus_poly) == rem((x + 1)^10, modulus_poly)
        @test PowerMod(x + 1, -1, modulus_poly) == invg

        @test truncate(f, 3) == poly(F, coefficients(f)[1:3])
        @test MulTrunc(f, g, 4) == truncate(f * g, 4)
        h = x^3 + 2x + 1
        invh = InvTrunc(h, 8)
        @test truncate(h * invh, 8) == one(h)
        @test degree(random_poly(F, 6; rng = MersenneTwister(1))) <= 5

        roots = F.([0, 1, -1])
        fromroots = BuildFromRoots(F, roots)
        @test all(iszero(evaluate(fromroots, r)) for r in roots)
        xs = F.([0, 1])
        ys = F.([1, 3])
        interpolated = interpolate(F, xs, ys)
        @test evaluate(interpolated, xs) == ys
        @test_throws ArgumentError interpolate(F, [0, 0], [1, 2])

        mutable_poly = x + 1
        setcoeff!(mutable_poly, 4, 2)
        @test coeff(mutable_poly, 4) == 2
        setcoeff!(mutable_poly, 4, 0)
        @test degree(mutable_poly) == 1
    end

    # ZZX-style construction and exact coefficient arithmetic reuse the same dense core.
    f = poly(BigInt, [1, 2, 3])
    g = poly(BigInt, [2, -1])
    @test coefficients(f * g) == BigInt[2, 3, 4, -3]
    @test evaluate(f, BigInt(2)) == 17
end

@testset "Polynomial factorization and irreducibility" begin
    rng = MersenneTwister(2026)
    for p in (2, 3, 5, 7, 11)
        F = zzp(p)
        x = xpoly(F)
        f = (x + 1)^2 * (x + 2) * (x^2 + x + 1)
        fac = factor(f; rng = rng)
        @test reconstruct(fac) == f
        @test all(is_irreducible(g) for (g, _) in fac.factors)
        @test all(m >= 1 for (_, m) in fac.factors)
        sf = SquareFreeDecomp(f)
        @test prod(g^m for (g, m) in sf; init = poly(F, 1)) == monic(f)
        roots = FindRoots((x + 1)^2 * (x - 2); rng = rng)
        @test length(roots) == 3
        @test count(r -> r == F(-1), roots) == (p == 3 ? 3 : 2)
        @test FindRoot(x + 3; rng = rng) == F(-3)
        @test FindRoot(x^2 + 1; rng = rng) === nothing || p % 4 == 1 || p == 2
        for d in 1:4
            irreducible = BuildIrred(F, d; rng = rng)
            @test degree(irreducible) == d
            @test is_irreducible(irreducible)
        end
        @test !is_irreducible((x + 1)^2)
        constant = factor(poly(F, -1); rng = rng)
        @test reconstruct(constant) == poly(F, -1)
    end

    # A nontrivial p-th-power branch of square-free factorization.
    F = zzp(3)
    x = xpoly(F)
    g = x^2 + 1
    f = g^3 * (x + 1)^2
    fac = factor(f; rng = rng)
    @test reconstruct(fac) == f
    multiplicities = Dict(coefficients(h) => m for (h, m) in fac.factors)
    @test sort([m for (_, m) in fac.factors]) == [2, 3]

    # The same algorithms work with arbitrary-precision prime moduli.
    FB = ZZP(big(2)^61 - 1)
    xb = xpoly(FB)
    fb = (xb + 1)^2 * (xb^2 + xb + 1)
    facb = factor(fb; rng = rng)
    @test reconstruct(facb) == fb
    @test_throws ArgumentError factor(poly(ZZP(12), [1, 1]))
end
