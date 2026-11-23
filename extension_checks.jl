using Test, Random, NTL
@testset "generic vs extension-field factorization" begin
    for (p,n) in ((2,3),(3,2),(5,2))
        E,a=conway_field(p,n)
        x=xpoly(E)
        h=BuildIrred(E,4;rng=MersenneTwister(729))
        for f in (h,h*(x+a)^2,h^3*(x+a)^2)
            z=p==2 ? GF2EXFactoring(f;rng=MersenneTwister(12)) : ZZ_pEXFactoring(f;rng=MersenneTwister(12))
            generic=factor(f;rng=MersenneTwister(82))
            @test reconstruct(z)==reconstruct(generic)==f
            @test length(generic.factors)==length(z.factors)
            @test all(t in generic.factors for t in z.factors)
        end
    end
    p=big(2)^127-1
    F=BigModRing(p)
    u=xpoly(F)
    E=extension_field(F,u^2+1)
    x=xpoly(E)
    a=gen(E)
    f=(x+a)^2
    @test (x+a,2) in ZZ_pEXFactoring(f).factors
end
