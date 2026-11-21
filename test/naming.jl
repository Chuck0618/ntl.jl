@testset "Public naming and compatibility aliases" begin
    @test ModRing === ZZPField
    @test ModRing === ZZP
    @test ModInt === ZZp
    @test ModInt === ZZ_p
    @test ExtensionElem === ZZpE

    F = ModRing(97)
    a = ModInt(F, 3)
    @test a == F(3)
    @test a isa ModInt
    @test ZZP(97) == F
    @test ZZPField(97) == F
    @test ZZp(F, 3) == a
    @test ZZ_p(F, 3) == a

    base = zzp(5)
    E = extension_field(base, xpoly(base)^2 + 2)
    b = ExtensionElem(E, [1, 1])
    @test b == E([1, 1])
    @test b isa ExtensionElem
    @test ZZpE(E, [1, 1]) == b
end
