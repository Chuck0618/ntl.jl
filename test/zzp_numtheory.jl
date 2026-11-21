@testset "Arbitrary-precision modular integers" begin
    for n in (2, 7, 12, 97)
        F = ModRing(n)
        @test modulus(F) == n
        @test characteristic(F) == n
        for a in 0:(n - 1), b in 0:(n - 1)
            x, y = F(a), F(b)
            @test rep(x + y) == mod(a + b, n)
            @test rep(x - y) == mod(a - b, n)
            @test rep(x * y) == mod(a * b, n)
            @test rep(-x) == mod(-a, n)
            if gcd(b, n) == 1
                @test y * inv(y) == one(F)
                @test x / y == x * inv(y)
            end
        end
        @test zero(F) == F(0)
        @test one(F) == F(1)
        @test oneunit(F) == F(1)
        @test F(-1) == F(n - 1)
        @test parent(F(2)) == F
        @test BigInt(F(2)) == mod(2, n)
    end
    F = ModRing(big(10)^60 + 37)
    a = F(-1)
    @test a * a == one(F)
    @test a^2 == one(F)
    @test a^(-1) == a
    @test rep(a) == modulus(F) - 1
    @test_throws ArgumentError F(1) + ModRing(7)(1)
    @test_throws DomainError inv(ModRing(12)(6))
    @test ModRing(12)(5)^-1 == ModRing(12)(5)
    @test VectorConv(ModRing(7), [-1, 8]) == ModRing(7).([-1, 8])
    r1, r2 = MersenneTwister(3), MersenneTwister(3)
    @test rand(r1, F, 10) == rand(r2, F, 10)
    @test all(x -> 0 <= rep(x) < modulus(F), rand(MersenneTwister(4), F, 20))
    @test random_zzp(F; rng = MersenneTwister(5)) isa ModInt
    @test length(VectorRandom(F, 5; rng = MersenneTwister(6))) == 5
end

@testset "Primality, random integers, and integer factorization" begin
    @test ProbPrime(big(2)^127 - 1)
    @test !ProbPrime(3 * (big(2)^127 - 1))
    for n in (4, 9, 25, 341, 561, 645, 1105, 1729)
        @test !ProbPrime(n)
    end
    @test MillerWitness(341, 2)
    @test !MillerWitness(17, 2)
    @test NextPrime(11) == 11
    @test NextPrime(12) == 13
    @test PrevPrime(13) == 13
    @test PrevPrime(12) == 11
    @test PrevPrime(2) == 2
    @test_throws DomainError PrevPrime(1)
    rng = MersenneTwister(2026)
    for bits in (1, 2, 8, 65, 130)
        x = RandomLen(bits; rng = rng)
        @test ndigits(x, base = 2) == bits
        @test 0 <= RandomBits(bits; rng = rng) < (big(1) << bits)
    end
    p = GenPrime(32; rng = rng)
    @test ndigits(p, base = 2) == 32
    @test ProbPrime(p)
    @test 0 <= RandomBnd(100; rng = rng) < 100
    @test_throws DomainError RandomBnd(0)

    for n in (1, 2, 12, 97, 99991 * 10007, -360)
        fac = FactorInt(n)
        value = prod(p^e for (p, e) in fac; init = BigInt(1))
        @test value == n
        @test all(ProbPrime(abs(p)) || p == -1 for (p, _) in fac)
    end
    @test EulerPhi(1) == 1
    @test EulerPhi(36) == 12
    @test EulerPhi(97) == 96
    @test_throws DomainError FactorInt(0)
end

@testset "CRT, modular square roots, discrete logarithms, and primitive roots" begin
    x, M = CRT([2, 3, 2], [3, 5, 7])
    @test (x, M) == (23, 105)
    @test all(x % m == r for (r, m) in zip([2, 3, 2], [3, 5, 7]))
    @test CRT([5], [1]) == (0, 1)
    @test_throws ArgumentError CRT([1, 2], [4, 6])
    @test_throws DimensionMismatch CRT([1], [2, 3])

    for p in (2, 3, 5, 7, 13, 17, 41, 97)
        for a in 0:(p - 1)
            root = SqrRootMod(a, p)
            brute = findfirst(y -> y^2 % p == a, 0:(p - 1))
            if brute === nothing
                @test root === nothing
            else
                @test root !== nothing
                @test mod(root^2, p) == a
            end
        end
    end
    @test_throws ArgumentError SqrRootMod(2, 15)

    @test DiscreteLog(2, 8, 13) == 3
    @test DiscreteLog(2, 1, 13) == 0
    @test DiscreteLog(2, 0, 13) === nothing
    g = PrimitiveRoot(97)
    @test g == 5
    seen = Set(powermod(g, e, 97) for e in 0:95)
    @test length(seen) == 96
    for e in (0, 1, 2, 17, 95)
        x = DiscreteLog(g, powermod(g, e, 97), 97)
        @test x !== nothing
        @test powermod(g, x, 97) == powermod(g, e, 97)
    end
end
