@testset "Public naming and compatibility aliases" begin
    @test PrimeField === zzp
    @test BigModRing === ModRing
    @test BigModRing === ZZPField
    @test BigModRing === ZZP
    @test ModInt === ZZp
    @test ModInt === ZZ_p
    @test ExtensionElem === ZZpE

    F = BigModRing(97)
    a = ModInt(F, 3)
    @test a == F(3)
    @test a isa ModInt
    @test ModRing(97) == F
    @test ZZP(97) == F
    @test ZZPField(97) == F
    @test ZZp(F, 3) == a
    @test ZZ_p(F, 3) == a

    base = PrimeField(5)
    E = extension_field(base, xpoly(base)^2 + 2)
    b = ExtensionElem(E, [1, 1])
    @test b == E([1, 1])
    @test b isa ExtensionElem
    @test ZZpE(E, [1, 1]) == b

    P = PrimeField(5)
    @test P === PrimeField{5}
    @test zzp(5) === P
    @test zzp{5}(3) == P(3)
    @test P(3) isa PrimeField
end
