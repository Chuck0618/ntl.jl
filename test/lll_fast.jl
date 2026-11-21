@testset "Fast LLL with exact verification" begin
    A = BigInt[1 1 1; -1 0 2; 3 5 6]
    fast = lll(A; mode = :fast)
    exact = lll(A; mode = :exact)

    for result in (fast, exact)
        B, U = result.basis, result.transform
        @test U * A == B
        @test abs(det_exact(U)) == 1
        @test result.verified
        mu, norms = gram_schmidt_exact(B)
        for i in 1:size(B, 1), j in 1:(i - 1)
            @test abs(mu[i, j]) <= BigInt(1) // BigInt(2)
        end
        for k in 2:size(B, 1)
            @test norms[k] >= (BigInt(3) // BigInt(4) - mu[k, k - 1]^2) * norms[k - 1]
        end
    end
    @test fast.mode == :fast
    @test exact.mode == :exact

    unchecked = lll(A; mode = :fast, verify = false)
    @test !unchecked.verified
    @test unchecked.transform * A == unchecked.basis
    @test_throws ArgumentError lll(A; mode = :unknown)

    B = BigInt[4 1 0; 1 5 2; 0 2 6]
    result = lll(B; mode = :fast)
    @test result.transform * B == result.basis
    @test result.verified
end
