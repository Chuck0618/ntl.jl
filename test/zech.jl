@testset "Zech extension arithmetic" begin
    for p in (2, 3, 5, 7)
        F = PrimeField(p)
        m = p == 2 ? xpoly(F)^2 + xpoly(F) + 1 : (p == 3 ? xpoly(F)^2 + xpoly(F) + 2 : (p == 5 ? xpoly(F)^2 + 2 : xpoly(F)^2 + 1))
        E = extension_field(F, m)
        Z = zech_field(E)
        q = cardinality(Z)
        @test q == p^2
        @test length(Z) == q
        @test length(collect(Z)) == q
        for a in collect(Z)
            @test a + zero(Z) == a
            @test a * one(Z) == a
            @test a - a == zero(Z)
            if !iszero(a)
                @test a * inv(a) == one(Z)
                @test a^(q - 1) == one(Z)
                @test a^q == a
            end
        end
        a = rand(MersenneTwister(p), Z)
        b = rand(MersenneTwister(p + 1), Z)
        @test muladd(a, b, b) == a * b + b
        @test isbitstype(typeof(a))
        @test extension_element(a * b) == extension_element(a) * extension_element(b)
        @test extension_element(a + b) == extension_element(a) + extension_element(b)
        @test extension_element(a - b) == extension_element(a) - extension_element(b)
        if !iszero(a)
            @test extension_element(inv(a)) == inv(extension_element(a))
            @test extension_element(a^3) == extension_element(a)^3
        end
        @test rand(MersenneTwister(1), Z, 5) == rand(MersenneTwister(1), Z, 5)
        @test Int(a) == rep(a)
        @test_throws DomainError inv(zero(Z))
    end

    F = PrimeField(3)
    E = extension_field(F, xpoly(F)^2 + xpoly(F) + 2)
    @test_throws ArgumentError zech_field(E; primitive = zero(E))
    Fbig = PrimeField(17)
    Ebig = extension_field(Fbig, BuildIrred(Fbig, 4; rng = MersenneTwister(2)))
    @test_throws ArgumentError zech_field(Ebig)
end
