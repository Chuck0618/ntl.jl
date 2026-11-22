@testset "Canonical public names and legacy cleanup" begin
    P = PrimeField(5)
    @test P === PrimeField{5}
    @test P(3) isa PrimeField

    F = BigModRing(97)
    a = ModInt(F, 3)
    @test a == F(3)
    @test a isa ModInt

    E = extension_field(P, xpoly(P)^2 + 2)
    b = ExtensionElem(E, [1, 1])
    @test b == E([1, 1])
    @test b isa ExtensionElem

    for name in (:zzp, :ModRing, :ZZPField, :ZZP, :ZZp, :ZZ_p, :ZZpE,
                 :zzp_Min, :zzp_Max, :random_zzp, :LLL)
        @test !isdefined(NTL, name)
    end
end
