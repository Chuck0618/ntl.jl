using Test, Random, NTL
rng = MersenneTwister(423)
d = big(3)//4
reference(A,B,U,d) = NTL._lll_verify(A,big.(B),big.(U),d)
verify_reduction(A,B,U,d) = NTL._lll_verify_reduction(A,B,U,d,big.(B),big.(U))
@testset "exact certificates" begin
    for n in 2:14, extra in (0,3), trial in 1:4
        A = rand(rng,-8:8,n,n+extra)
        r = lll(A)
        @test eltype(r.basis) == eltype(r.transform) == BigInt
        @test r.verified && reference(A,r.basis,r.transform,d)
        B,U = Int.(r.basis),Int.(r.transform)
        for q in (d, big(99)//100)
            @test verify_reduction(A,B,U,q) == reference(A,B,U,q)
        end
        B[1,1] += 1
        @test !verify_reduction(A,B,U,d)
    end
    I2 = [1 0;0 1]
    for (B,U) in (([2 0;0 2],[2 0;0 2]), ([1 0;2 1],I2), ([3 0;0 1],I2), ([1 0;1 0],I2))
        A = U == I2 ? B : I2
        @test !verify_reduction(A,B,U,d)
    end
    @test verify_reduction([0 1;1 0],I2,[0 1;1 0],d)
    @test verify_reduction([2 0;1 2],[2 0;1 2],I2,d)
    for scale in (1<<15,1<<21)   # wide Gram minors, and entry-size fallback
        B = [i==j ? scale : 0 for i=1:12,j=1:12]
        U = [Int(i==j) for i=1:12,j=1:12]
        @test verify_reduction(B,B,U,d)
    end
    @test verify_reduction(I2,I2,I2,((big(1)<<160)-1)//(big(1)<<160))
    wide = big.(I2) .* (big(1)<<75)
    @test lll(wide).verified && reference(wide,lll(wide).basis,lll(wide).transform,d)
    @test verify_reduction([2 0;0 1],[2 0;0 1],I2,d) == false
    # Exact Lovasz equality and an arbitrarily close strict failure.
    edge = [1 1;1 0]
    @test verify_reduction(edge,edge,I2,big(1)//2)
    @test !verify_reduction(edge,edge,I2,((big(1)<<60)+1)//(big(1)<<61))
    # Unimodular block matrix with large leading minors: Int Bareiss overflows.
    eye = [Int(i==j) for i in 1:5,j in 1:5]
    z = zeros(Int,5,5)
    Uwide = [100000eye eye; eye z]
    Awide = [z eye; eye -100000eye]
    @test verify_reduction(Awide,Awide*Uwide,Uwide,d)
    if !isempty(ARGS)
        values = parse.(Int,readlines(ARGS[1]))
        B = permutedims(reshape(values[1:144],12,12))
        U = permutedims(reshape(values[145:end],12,12))
        A = [i==j ? i+1 : j<=i ? i*(j+1)%17 : 0 for i=1:12,j=1:12]
        @test reference(A,B,U,d)
        @test verify_reduction(A,B,U,d)
    end

    # A alone may exceed Int even when B and U are small. Conversion must not
    # throw before the generic verifier gets control.
    Awide = (big(1)<<75) * BigInt[1 0; 0 1]
    Ibig = BigInt[1 0; 0 1]
    Bsmall, Usmall = Int[1 0; 0 1], Int[1 0; 0 1]
    @test verify_reduction(Awide, Awide, Ibig, d)
    @test !verify_reduction(Awide, Bsmall, Usmall, d)
end
