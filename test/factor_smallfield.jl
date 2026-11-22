@testset "Small-modulus arithmetic and factorization" begin
    rng = MersenneTwister(201)

    for p in (2, 3, 97, 157, 199, 251, 36587, 44203, 65521, 65537,
              8423879, 1_000_000_007, 2147483647, 4112838589)
        p > typemax(Int) && continue
        F = PrimeField(Int(p))
        if p <= 251
            for a in 0:(p - 1), b in 0:(p - 1)
                @test rep(F(a) * F(b)) == mod(big(a) * b, p)
                @test rep(F(a) + F(b)) == mod(big(a) + b, p)
            end
        else
            for _ in 1:2000
                a, b = rand(rng, 0:(p - 1)), rand(rng, 0:(p - 1))
                @test rep(F(a) * F(b)) == mod(big(a) * b, p)
                @test rep(F(a) + F(b)) == mod(big(a) + b, p)
            end
        end
    end

    for p in (2, 3, 5, 17, 101, 127, 131, 1009)
        F = PrimeField(p)
        x = xpoly(F)
        for d in (1, 2, 3, 4, 7, 8, 16, 17), trial in 1:4
            f = poly(F, vcat(rand(rng, 0:(p - 1), d), rand(rng, 1:(p - 1))))
            before = coefficients(f)
            result = NTL.factor(f; rng=rng)
            @test reconstruct(result) == f
            @test all(pair -> isone(LeadCoeff(first(pair))) &&
                              is_irreducible(first(pair)) && last(pair) > 0,
                      result.factors)
            @test all(i -> all(j -> result.factors[i][1] != result.factors[j][1],
                               1:(i - 1)), eachindex(result.factors))
            @test coefficients(f) == before
        end
        @test_throws DomainError NTL.factor(zero(x))
        @test isempty(NTL.factor(poly(F, 1)).factors)
        @test NTL.factor((x + 1)^8).factors == [(x + 1, 8)]
    end

    # Exhaust every monic quartic over F_2 and F_3. In particular, this covers
    # root-free quartics handled by _quartic_pair and compares against CanZass.
    for p in (2, 3)
        F = PrimeField(p)
        for A in Iterators.product(ntuple(_ -> 0:(p - 1), 4)...)
            f = poly(F, vcat(collect(A), [one(F)]))
            result = NTL.factor(f; rng=MersenneTwister(2026))
            reference = NTL._factor_generic(f, MersenneTwister(2026))
            @test reconstruct(result) == f
            @test Set((coefficients(g), m) for (g, m) in result.factors) ==
                  Set((coefficients(g), m) for (g, m) in reference.factors)
        end
    end

    for p in (31, 101)
        F = PrimeField(p)
        x = xpoly(F)
        f = (x + 1) * (x + 2) * (x^2 + 2)
        result = NTL.factor(f; rng=MersenneTwister(2026))
        reference = NTL._factor_generic(f, MersenneTwister(2026))
        @test reconstruct(result) == f
        @test Set((coefficients(g), m) for (g, m) in result.factors) ==
              Set((coefficients(g), m) for (g, m) in reference.factors)
    end
end
