using NTL
using Random
using Statistics
using Printf

function median_ms(f; reps = 8)
    f()
    times = [@elapsed f() for _ in 1:reps]
    1000 * Statistics.median(times)
end

F = PrimeField(998244353)
n = 600
f = random_poly(F, n; rng = MersenneTwister(2026))
g = random_poly(F, n; rng = MersenneTwister(2027))

reference = poly_mul(f, g; method = :schoolbook)
for method in (:schoolbook, :karatsuba, :ntt, :auto)
    result = poly_mul(f, g; method = method)
    @assert result == reference
    ms = median_ms(() -> poly_mul(f, g; method = method))
    @printf("%-12s %10.4f ms\n", method, ms)
end

# The default is the original schoolbook path. Enable degree-based acceleration
# only when desired; unsupported NTT moduli fall back to Karatsuba.
set_polynomial_acceleration!(true)
@assert polynomial_acceleration()
@assert poly_mul(f, g) == reference
set_polynomial_acceleration!(false)
@assert !polynomial_acceleration()

println("Polynomial acceleration example passed.")
