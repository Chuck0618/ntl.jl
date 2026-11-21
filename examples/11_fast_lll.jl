using NTL
using Statistics
using Printf

A = BigInt[1 1 1; -1 0 2; 3 5 6]
result = lll(A; mode = :fast, verify = true)
B, U = result.basis, result.transform

@assert U * A == B
@assert abs(det_exact(U)) == 1
@assert result.verified

println("mode = ", result.mode)
println("LLL-reduced basis = ", B)
println("unimodular transform = ", U)

for mode in (:fast, :exact)
    started = time()
    r = lll(A; mode = mode, verify = true)
    elapsed = 1000 * (time() - started)
    @printf("%-8s %10.4f ms (verified=%s)\n", mode, elapsed, r.verified)
    @assert r.transform * A == r.basis
end
