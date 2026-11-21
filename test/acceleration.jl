@testset "Polynomial multiplication acceleration" begin
    @test polynomial_acceleration()
    set_polynomial_acceleration!(false)
    @test !polynomial_acceleration()

    for p in (2, 3, 5, 101), n in (1, 5, 32, 80, 160)
        F = zzp(p)
        f = random_poly(F, n; rng = MersenneTwister(n + p))
        g = random_poly(F, n; rng = MersenneTwister(2n + p))
        reference = poly_mul(f, g; method = :schoolbook)
        @test poly_mul(f, g; method = :karatsuba) == reference
        @test poly_mul(f, g) == reference
        @test poly_mul(f, g; method = :auto) == reference
        set_polynomial_acceleration!(true)
        @test poly_mul(f, g) == reference
        set_polynomial_acceleration!(false)
    end

    F = zzp(97)
    f = random_poly(F, 320; rng = MersenneTwister(1))
    g = random_poly(F, 320; rng = MersenneTwister(2))
    reference = poly_mul(f, g; method = :schoolbook)
    @test poly_mul(f, g; method = :karatsuba) == reference
    @test poly_mul(f, g; method = :auto) == reference
    @test_throws ArgumentError poly_mul(f, g; method = :ntt)
    @test_throws ArgumentError poly_mul(f, g; method = :unknown)

    F = zzp(998244353)
    for n in (64, 300, 600)
        f = random_poly(F, n; rng = MersenneTwister(n + 1))
        g = random_poly(F, n; rng = MersenneTwister(n + 2))
        reference = poly_mul(f, g; method = :schoolbook)
        @test poly_mul(f, g; method = :ntt) == reference
        @test poly_mul(f, g; method = :auto) == reference
        set_polynomial_acceleration!(true)
        @test poly_mul(f, g) == reference
        set_polynomial_acceleration!(false)
    end
    set_polynomial_acceleration!(false)
end

@testset "Generalized FFT convolution and CRT reconstruction" begin
    for p in (2, 3, 5, 101, 997, 9973, 1000003)
        F = zzp(p)
        n = p == 1000003 ? 300 : 80
        a = rand(MersenneTwister(p + 1), F, n)
        b = rand(MersenneTwister(p + 2), F, n)
        @test fft_convolve(a, b) == convolve(a, b)
        f = poly(F, a)
        g = poly(F, b)
        reference = poly_mul(f, g; method = :schoolbook)
        @test poly_mul(f, g; method = :fft) == reference
        @test poly_mul(f, g; method = :auto) == reference
    end

    F = zzp(1000003)
    a = fill(F(1000002), 1000)
    b = fill(F(1000002), 1000)
    moduli = NTL._fft_prime_plan(NTL._fft_bound(a, b), length(a) + length(b) - 1)
    @test length(moduli) >= 2
    @test fft_convolve(a, b) == convolve(a, b)
    @test NTL._crt_coefficients([[2, 2], [3, 2], [2, 5]], [3, 5, 7], 2) == BigInt[23, 47]

    G = ZZP(12)
    a = G.([1, 5, 7, 11])
    b = G.([3, 2, 9])
    @test fft_convolve(a, b) == convolve(a, b)

    a = GF2.([1, 0, 1, 1])
    b = GF2.([1, 1, 0])
    @test fft_convolve(a, b) == convolve(a, b)
end
