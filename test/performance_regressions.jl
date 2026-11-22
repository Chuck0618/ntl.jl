@testset "Optimized field, polynomial, and lattice kernels" begin
    rng = MersenneTwister(1771)
    small_primes = (2, 101, 998244353, PRIME_FIELD_MAX)
    word_primes = Sys.WORD_SIZE == 64 ? (small_primes..., 3221225473, 2305843009213693951) : small_primes
    for p in word_primes
        F = PrimeField(p)
        for _ in 1:80
            a, b, c = rand(rng, 0:(p - 1), 3)
            @test rep(F(a) * F(b)) == mod(BigInt(a) * b, p)
            @test rep(muladd(F(a), F(b), F(c))) == mod(BigInt(a) * b + c, p)
            @test rep(power(F(a), typemax(UInt))) == powermod(BigInt(a), BigInt(typemax(UInt)), BigInt(p))
        end
        @test rep(power(F(p - 1), typemin(Int))) == 1
    end
    for n in (big(1), big(18), (big(1) << 128) + 51)
        F, G = BigModRing(n), BigModRing(n)
        a, b = F(n - 1), G(n - 1)
        @test rep(a + b) == mod(2n - 2, n)
        @test rep(a - b) == 0
        @test rep(a * b) == mod((n - 1)^2, n)
        @test rep(a^-3) == mod(n - 1, n)
        @test rep(a) == n - 1 && rep(b) == n - 1
    end
    @test_throws DomainError power(BigModRing(12)(6), -1)
    @test_throws ArgumentError BigModRing(12)(1) + BigModRing(15)(1)

    transform_primes = Sys.WORD_SIZE == 64 ? (17, 998244353, 3221225473, 4179340454199820289) : (17, 998244353)
    for p in transform_primes
        F = PrimeField(p)
        for (a, b) in ((1, 1), (2, 3), (7, 10))
            A, B = rand(rng, F, a), rand(rng, F, b)
            @test ntt_convolve(A, B) == convolve(A, B)
            @test intt(ntt(vcat(A, fill(zero(F), 16 - a)))) == vcat(A, fill(zero(F), 16 - a))
        end
    end
    for k in 1:6
        moduli = collect(FFT_PRIMES[1:k])
        product = prod(BigInt.(moduli))
        samples = BigInt[0, 1, product ÷ 2, product - 1]
        residues = [[Int(x % q) for x in samples] for q in moduli]
        @test NTL._crt_coefficients(residues, moduli, length(samples)) == samples
    end
    for (p, n, m) in ((2, 128, 32), (101, 320, 257), (65537, 128, 160),
                      (998244353, 144, 130))
        F = PrimeField(p)
        a = poly(F, vcat(rand(rng, 0:(p - 1), n - 1), [p - 1]))
        b = poly(F, fill(p - 1, m))
        @test a * b == poly_mul(a, b; method=:schoolbook)
        @test b * a == a * b
    end

    for p in (2, 101, PRIME_FIELD_MAX), d in (1, 2, 7, 16, 17)
        F = PrimeField(p)
        f = poly(F, rand(rng, 0:(p - 1), d + 3))
        m = poly(F, vcat(rand(rng, 0:(p - 1), d), [p - 1]))
        for e in (0, 1, 31, 101, big(1) << 70)
            @test PowerMod(f, e, m) == NTL._polynomial_power_binary(rem(f, m), e, m)
        end
        @test PowerMod(one(f), typemin(Int), m) == one(f)
    end

    for p in (2, 101, 1_000_000_007, PRIME_FIELD_MAX)
        F = PrimeField(p)
        for n in (1, 7, 12, 40)
            A = rand(rng, F, n, n)
            @test det(A) == NTL._det_generic(A, n)
            A[end, :] .= zero(F)
            @test iszero(det(A))
        end
    end

    for n in (2, 5, 12), width in (n, n + 3)
        B = BigInt.(rand(rng, -100:100, n, width))
        for i in 1:n
            B[i, i] += 10_000
        end
        mu, norms = gram_schmidt_exact(B)
        L = mu + [BigInt(i == j) for i in 1:n, j in 1:n]
        for i in 1:n, j in 1:n
            @test sum(L[i, k] * norms[k] * L[j, k] for k in 1:n) ==
                  sum(B[i, k] * B[j, k] for k in 1:width)
        end
        reduced = lll(B; delta=99//100)
        @test reduced.verified
        @test reduced.transform * B == reduced.basis
        @test abs(det_exact(reduced.transform)) == 1
    end
end
