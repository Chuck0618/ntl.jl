@testset "Exact matrices over supported fields" begin
    parents = Any[PrimeField(7), BigModRing(11), GF2]
    F5 = PrimeField(5)
    E = extension_field(F5, xpoly(F5)^2 + 2)
    push!(parents, E)
    for F in parents
        A = F.([2 3 1; 4 1 5; 3 7 2])
        b = F.([1, 2, 3])
        d = det(A)
        if !iszero(d)
            x = solve(A, b)
            @test A * x == b
            invA = matrix_inverse(A)
            @test A * invA == F.([1 0 0; 0 1 0; 0 0 1])
            @test det(invA) == inv(d)
        else
            @test_throws ArgumentError solve(A, b)
        end
        singular = F.([1 2 3; 2 4 6; 0 1 1])
        @test iszero(det(singular))
        @test rank(singular) == 2
        ns = nullspace(singular)
        @test length(ns) == 1
        @test all(iszero, singular * ns[1])
        R, pivots = rref(singular)
        @test pivots == [1, 2]
        @test_throws ArgumentError matrix_inverse(singular)
    end
    @test det_exact([2 3; 4 5]) == -2
    @test det_exact(big.([1 2 3; 0 4 5; 0 0 6])) == 24
    @test det_exact([1 2; 2 4]) == 0
    @test det_exact(reshape([7], 1, 1)) == 7
end

@testset "Exact LLL lattice reduction" begin
    A = BigInt[1 1 1; -1 0 2; 3 5 6]
    result = lll(A)
    B, U = result.basis, result.transform
    @test U * A == B
    @test abs(det_exact(U)) == 1
    mu, norms = gram_schmidt_exact(B)
    for i in 1:size(B, 1), j in 1:(i - 1)
        @test abs(mu[i, j]) <= BigInt(1) // BigInt(2)
    end
    for k in 2:size(B, 1)
        @test norms[k] >= (BigInt(3) // BigInt(4) - mu[k, k - 1]^2) * norms[k - 1]
    end
    @test LLL(A).basis == B
    identity_result = lll([1 0; 0 1])
    @test identity_result.basis == BigInt[1 0; 0 1]
    @test identity_result.transform == BigInt[1 0; 0 1]
    @test lll(reshape(BigInt[7], 1, 1); delta = 0.9).basis == reshape(BigInt[7], 1, 1)
    @test_throws DomainError lll(A; delta = 0.2)
    @test_throws ArgumentError lll([1 2; 2 4])
end

@testset "NTT and convolution" begin
    for F in (PrimeField(17), BigModRing(17))
        g = primitive_root(F)
        @test g^(16) == one(F)
        @test g^8 != one(F)
        omega = root_of_unity(8, F)
        @test omega^8 == one(F)
        @test omega^4 != one(F)
        a = rand(MersenneTwister(22), F, 8)
        @test intt(ntt(a)) == a
        b = rand(MersenneTwister(23), F, 5)
        @test ntt_convolve(a, b) == convolve(a, b)
        @test ntt([F(3)]) == [F(3)]
        @test intt([F(3)]) == [F(3)]
    end
    F = PrimeField(97)
    a = rand(MersenneTwister(24), F, 13)
    b = rand(MersenneTwister(25), F, 9)
    @test ntt_convolve(a, b) == convolve(a, b)
    @test isempty(convolve(F[], F[]))
    @test_throws ArgumentError ntt(F.(1:3))
    @test_throws ArgumentError root_of_unity(5, PrimeField(17))
end
