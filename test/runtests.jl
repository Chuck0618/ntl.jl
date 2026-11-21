using Test
using Random
using NTL

ref_mod(a, n) = Int(mod(BigInt(a), BigInt(n)))

function ref_powermod(a, e, n)
    base = mod(BigInt(a), BigInt(n))
    exponent = BigInt(e)
    modulus_big = BigInt(n)
    if exponent < 0
        base = invmod(base, modulus_big)
        exponent = -exponent
    end
    result = mod(BigInt(1), modulus_big)
    while exponent > 0
        if isodd(exponent)
            result = mod(result * base, modulus_big)
        end
        exponent >>= 1
        if exponent > 0
            base = mod(base * base, modulus_big)
        end
    end
    Int(result)
end

@testset "NTL single-word prime fields" begin
    @testset "modulus discovery and validation" begin
        @test zzp_Min == 2
        @test NTL._is_prime(zzp_Max)
        @test zzp_Max <= typemax(Int)
        @test all(!NTL._is_prime(n) for n in (zzp_Max + 1):typemax(Int))
        @test zzp(7) === zzp{7}
        @test zzp(big(97)) === zzp{97}
        @test zzp(zzp_Max) === zzp{zzp_Max}
        for p in (2, 3, 5, 7, 11, 37, 41, 97, 65537)
            @test zzp(p) === zzp{p}
        end
        for n in (-7, -1, 0, 1)
            @test_throws DomainError zzp(n)
        end
        for n in (4, 6, 9, 21, 25, 341, 645, 1105)
            @test_throws ArgumentError zzp(n)
        end
        @test_throws DomainError zzp(BigInt(typemax(Int)) + 1)
        @test_throws ArgumentError zzp{4}(1)
    end

    @testset "construction, conversion, and metadata" begin
        F = zzp(7)
        @test modulus(F) == 7
        @test modulus(F(3)) == 7
        @test characteristic(F) == 7
        @test characteristic(F(3)) == 7
        @test rep(F(0)) == 0
        @test rep(F(11)) == 4
        @test rep(F(-1)) == 6
        @test rep(F(-15)) == 6
        @test rep(F(big(10)^100 + 3)) == ref_mod(big(10)^100 + 3, 7)
        @test rep(zzp{7}()) == 0
        @test F(F(4)) === F(4)
        @test Int(F(6)) == 6
        @test convert(BigInt, F(6)) == 6
        @test convert(F, -9) == F(5)
        @test convert(F, F(5)) === F(5)
        @test_throws ArgumentError convert(F, zzp{11}(5))
        @test_throws ArgumentError F(zzp{11}(5))
        @test zero(F) == F(0)
        @test zero(F(3)) == F(0)
        @test one(F) == F(1)
        @test one(F(3)) == F(1)
        @test oneunit(F) == F(1)
        @test oneunit(F(3)) == F(1)
        @test iszero(F(14))
        @test !iszero(F(3))
        @test isone(F(8))
        @test !isone(F(0))
        @test info(zzp).minimum_modulus == 2
        @test info(zzp).maximum_modulus == zzp_Max
        @test info(F).modulus == 7
        @test info(F).storage == Int
        @test info(F(6)).representative == 6
        @test info(F(6)).field === F
        @test sprint(show, F(3)) == "zzp{7}(3)"
        @test sprint(show, F) == "zzp{7}"
        x, y = promote(F(2), 3)
        @test x isa F
        @test y isa F
        @test y == F(3)
    end

    @testset "addition, subtraction, and negation against BigInt" begin
        for p in (2, 3, 5, 7, 17, 101)
            F = zzp(p)
            for a in 0:(p - 1), b in 0:(p - 1)
                x, y = F(a), F(b)
                @test rep(x + y) == ref_mod(a + b, p)
                @test rep(x - y) == ref_mod(a - b, p)
                @test rep(add(x, y)) == ref_mod(a + b, p)
                @test rep(sub(x, y)) == ref_mod(a - b, p)
                @test rep(x + b) == ref_mod(a + b, p)
                @test rep(a + y) == ref_mod(a + b, p)
                @test rep(x - b) == ref_mod(a - b, p)
                @test rep(a - y) == ref_mod(a - b, p)
            end
            for a in 0:(p - 1)
                @test rep(-F(a)) == ref_mod(-a, p)
                @test rep(negate(F(a))) == ref_mod(-a, p)
                @test F(a) + (-F(a)) == zero(F)
            end
        end
        F = zzp(7)
        @test F(3) + big(9) == F(5)
        @test big(-8) - F(2) == F(4)
        @test_throws ArgumentError F(1) + zzp{11}(2)
        @test_throws ArgumentError F(1) - zzp{11}(2)
    end

    @testset "multiplication, squaring, and division against BigInt" begin
        for p in (2, 3, 5, 7, 17, 101)
            F = zzp(p)
            for a in 0:(p - 1), b in 0:(p - 1)
                x, y = F(a), F(b)
                @test rep(x * y) == ref_mod(BigInt(a) * b, p)
                @test rep(mul(x, y)) == ref_mod(BigInt(a) * b, p)
                @test rep(x * b) == ref_mod(BigInt(a) * b, p)
                @test rep(a * y) == ref_mod(BigInt(a) * b, p)
                if b != 0
                    @test rep(x / y) == ref_mod(BigInt(a) * invmod(BigInt(b), BigInt(p)), p)
                    @test rep(divide(x, y)) == ref_mod(BigInt(a) * invmod(BigInt(b), BigInt(p)), p)
                    @test y * inv(y) == one(F)
                end
            end
            for a in 0:(p - 1)
                @test rep(sqr(F(a))) == ref_mod(BigInt(a)^2, p)
            end
        end
        F = zzp(7)
        @test F(3) / 2 == F(5)
        @test 3 / F(2) == F(5)
        @test F(3) * big(8) == F(3)
        @test_throws DomainError inv(F(0))
        @test_throws DomainError F(1) / F(0)
        @test_throws DomainError F(1) / 7
        @test_throws ArgumentError F(1) * zzp{11}(1)
        @test_throws ArgumentError F(1) / zzp{11}(1)
    end

    @testset "low-level modular functions" begin
        for n in (1, 2, 3, 7, 12, 97, 1009)
            for a in 0:(n - 1), b in 0:(n - 1)
                @test AddMod(a, b, n) == ref_mod(a + b, n)
                @test SubMod(a, b, n) == ref_mod(a - b, n)
                @test MulMod(a, b, n) == ref_mod(BigInt(a) * b, n)
                @test SqrMod(a, n) == ref_mod(BigInt(a)^2, n)
            end
            for a in 0:(n - 1)
                @test NegateMod(a, n) == ref_mod(-a, n)
            end
        end
        @test InvMod(3, 7) == 5
        @test InvMod(10, 7) == 5
        @test InvMod(-4, 7) == 5
        @test InvMod(typemin(Int), 7) == ref_mod(invmod(BigInt(typemin(Int)), BigInt(7)), 7)
        @test_throws DomainError InvMod(3, 12)
        @test InvModStatus(5, 12) == (true, 5)
        @test InvModStatus(6, 12) == (false, 0)
        @test InvModStatus(0, 7) == (false, 0)
        @test_throws DomainError InvMod(6, 12)
        @test_throws DomainError InvMod(0, 7)
        @test_throws DomainError InvMod(1, 1)
        @test_throws DomainError AddMod(0, 0, 0)
        @test_throws DomainError AddMod(-1, 0, 7)
        @test_throws DomainError SubMod(0, 7, 7)
        @test_throws DomainError NegateMod(7, 7)
        @test_throws DomainError MulMod(0, -1, 7)
    end

    @testset "power and extreme exponents" begin
        for p in (2, 3, 5, 7, 17, 101)
            for a in 0:(p - 1), e in 0:12
                @test PowerMod(a, e, p) == ref_powermod(a, e, p)
            end
            for a in 1:(p - 1), e in -8:-1
                @test PowerMod(a, e, p) == ref_powermod(a, e, p)
            end
            F = zzp(p)
            for a in 0:(p - 1), e in 0:8
                @test F(a)^e == F(ref_powermod(a, e, p))
                @test power(F(a), e) == F(ref_powermod(a, e, p))
            end
            @test PowerMod(2, big(10)^30 + 7, p) == ref_powermod(2, big(10)^30 + 7, p)
        end
        F = zzp(7)
        @test F(0)^0 == one(F)
        @test F(3)^-1 == inv(F(3))
        @test F(3)^big(-17) == F(ref_powermod(3, -17, 7))
        @test PowerMod(3, typemin(Int), 97) == ref_powermod(3, typemin(Int), 97)
        @test PowerMod(-3, 9, 7) == ref_powermod(-3, 9, 7)
        @test PowerMod(0, 0, 1) == 0
        @test PowerMod(0, 5, 1) == 0
        @test_throws DomainError PowerMod(0, -1, 7)
        @test_throws DomainError PowerMod(2, 1, 0)
    end

    @testset "largest supported modulus and overflow boundaries" begin
        p = zzp_Max
        F = zzp(p)
        samples = unique((0, 1, 2, p ÷ 2, p - 3, p - 2, p - 1))
        for a in samples, b in samples
            @test AddMod(a, b, p) == ref_mod(BigInt(a) + b, p)
            @test SubMod(a, b, p) == ref_mod(BigInt(a) - b, p)
            @test MulMod(a, b, p) == ref_mod(BigInt(a) * b, p)
            @test rep(F(a) + F(b)) == ref_mod(BigInt(a) + b, p)
            @test rep(F(a) - F(b)) == ref_mod(BigInt(a) - b, p)
            @test rep(F(a) * F(b)) == ref_mod(BigInt(a) * b, p)
        end
        @test MulMod(p - 1, p - 1, p) == 1
        @test AddMod(p - 1, 1, p) == 0
        @test SubMod(0, 1, p) == p - 1
        @test PowerMod(p - 1, typemax(Int), p) == ref_powermod(p - 1, typemax(Int), p)
    end

    @testset "equality and field separation" begin
        F, G = zzp(7), zzp(11)
        @test F(3) == F(10)
        @test F(3) != F(4)
        @test F(3) == 3
        @test F(3) == 10
        @test F(3) == -4
        @test 10 == F(3)
        @test big(-4) == F(3)
        @test isequal(F(3), F(10))
        @test isequal(F(3), 10)
        @test isequal(10, F(3))
        @test F(3) != G(3)
        @test !isequal(F(3), G(3))
        @test hash(F(3)) == hash(F(10))
        @test swap(F(1), F(2)) == (F(2), F(1))
    end

    @testset "random elements and vector utilities" begin
        F = zzp(97)
        r1, r2 = MersenneTwister(2026), MersenneTwister(2026)
        a = rand(r1, F, 100)
        b = rand(r2, F, 100)
        @test a == b
        @test length(a) == 100
        @test all(x -> 0 <= rep(x) < 97, a)
        @test size(rand(MersenneTwister(1), F, 2, 3)) == (2, 3)
        @test random_zzp(F; rng = MersenneTwister(9)) isa F
        v = VectorRandom(F, 20; rng = MersenneTwister(10))
        @test length(v) == 20
        @test all(x -> x isa F, v)
        @test isempty(VectorRandom(F, 0))
        @test_throws DomainError VectorRandom(F, -1)

        converted = VectorConv(F, [-1, 0, 98, big(10)^30])
        @test converted == [F(-1), F(0), F(98), F(big(10)^30)]
        @test VectorMulMod([0, 1, 2, 96], 5, 97) == [0, 5, 10, 92]
        out = zeros(Int, 4)
        @test VectorMulMod!(out, [0, 1, 2, 96], 5, 97) === out
        @test out == [0, 5, 10, 92]
        aliased = [1, 2, 3]
        @test VectorMulMod!(aliased, aliased, 2, 7) === aliased
        @test aliased == [2, 4, 6]
        @test_throws DimensionMismatch VectorMulMod!(zeros(Int, 2), [1], 2, 7)
        @test_throws DomainError VectorMulMod([7], 2, 7)
    end
end

include("zzp_numtheory.jl")
include("polynomial.jl")
include("binary_extension.jl")
include("matrix_lattice_transform.jl")
include("acceleration.jl")
include("lll_fast.jl")
include("zech.jl")


@testset "zzp finite-field iteration" begin
    F = zzp(5)
    @test length(F) == 5
    @test eltype(F) === F
    @test Base.IteratorSize(F) == Base.HasLength()
    @test Base.IteratorEltype(F) == Base.HasEltype()
    @test first(F) == F(0)
    @test last(F) == F(4)
    @test collect(F) == F.(0:4)
    @test [rep(a) for a in F] == [0, 1, 2, 3, 4]
    count = 0
    total = zero(F)
    for a in F
        count += 1
        total += a
    end
    @test count == 5
    @test total == F(10)
    @test sum(F) == F(10)
    @test collect(zzp(2)) == zzp{2}.([0, 1])
    @test all(i - 1 == rep(a) for (i, a) in enumerate(F))
    @test_throws ArgumentError collect(zzp{4})
    @test_throws DomainError length(zzp{1})
end
