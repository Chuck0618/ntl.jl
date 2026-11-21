@testset "Public naming and compatibility aliases" begin
    @test ZZp === ZZ_p
    @test ZZpE === ExtensionElem

    F = ZZP(97)
    a = ZZp(F, 3)
    @test a == F(3)
    @test a isa ZZp
    @test a isa ZZ_p

    base = zzp(5)
    E = extension_field(base, xpoly(base)^2 + 2)
    b = ZZpE(E, [1, 1])
    @test b == E([1, 1])
    @test b isa ZZpE
    @test b isa ExtensionElem
end
