using NTL
using Random
using Statistics
using Printf

function median_ms(f; reps = 8)
    f()
    times = [@elapsed f() for _ in 1:reps]
    1000 * Statistics.median(times)
end

# This prime does not support the required transform length, so generalized FFT
# computes over compatible NTT primes and reconstructs exact coefficients.
F = zzp(101)
f = random_poly(F, 1000; rng = MersenneTwister(10))
g = random_poly(F, 1000; rng = MersenneTwister(20))
reference = poly_mul(f, g; method = :schoolbook)
@assert poly_mul(f, g; method = :fft) == reference
@assert poly_mul(f, g) == reference

for method in (:schoolbook, :karatsuba, :fft, :auto)
    ms = median_ms(() -> poly_mul(f, g; method = method))
    @printf("%-12s %10.4f ms\n", method, ms)
end

# A same-prime fast path avoids reconstruction when the original modulus already
# supports the transform.
G = zzp(998244353)
a = rand(MersenneTwister(1), G, 600)
b = rand(MersenneTwister(2), G, 600)
@assert fft_convolve(a, b) == convolve(a, b)

println("Generalized FFT example passed.")
