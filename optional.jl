using Test, NTL, Random
import Nemo, Polynomials, BitIntegers
rng=MersenneTwister(2094)
@testset "FLINT and optional companions" begin
    R,x=Nemo.polynomial_ring(Nemo.ZZ,:x)
    for m in 1:5,n in 1:5,k in 1:6
        A=rand(rng,-20:20,m,n)
        k==2 && (A[1,:].=0)
        a=Nemo.matrix(Nemo.ZZ,A)
        H=hnf(A;transform=false)
        @test Nemo.matrix(Nemo.ZZ,H)==Nemo.hnf(a)
        if m==n
            g=Nemo.charpoly(R,a)
            @test coefficients(charpoly(A))==[BigInt(Nemo.coeff(g,i)) for i in 0:n]
        end
    end
    for i in 1:12
        f=poly(BigInt,rand(rng,-9:9,6))
        g=poly(BigInt,rand(rng,-9:9,4))
        a=R(Nemo.ZZ.(coefficients(f)))
        b=R(Nemo.ZZ.(coefficients(g)))
        @test resultant(f,g)==BigInt(Nemo.resultant(a,b))
        @test discriminant(f)==BigInt(Nemo.discriminant(a))
        @test reconstruct(factor_integer_poly(f))==f
        z=Polynomials.Polynomial(f)
        @test ModPoly(z)==f
        @test ModPoly(z*Polynomials.Polynomial(g))==f*g
        @test z(7)==evaluate(f,7)
    end
    q=BigModRing(big(2)^195)
    v=packedfieldvec(q,[1,-1];storage=BitIntegers.Int256)
    @test v[2]==q(-1) && sizeof(v.data)==64
    @test_throws ArgumentError packedfieldvec(q,[1];storage=Int32)
    F=PrimeField(17)
    @test ModPoly(Polynomials.Polynomial(poly(F,[1,3])))==poly(F,[1,3])
end
@testset "external polynomial entry points and Conway source data" begin
    p=Polynomials.Polynomial(BigInt[1,3,2],:t)
    q=Polynomials.Polynomial(BigInt[1,1],:t)
    @test NTL.DivRem(p,q)==(Polynomials.Polynomial(BigInt[1,2],:t),zero(q))
    @test NTL.GCD(p,q)==q
    @test NTL.derivative(p)==Polynomials.derivative(p)
    @test NTL.evaluate(p,2)==p(2)
    @test reconstruct(NTL.factor(p))==ModPoly(p)
    @test NTL.discriminant(p)==1 && NTL.resultant(p,q)==0
    @test first(NTL.subresultants(p,q))==zero(p)
    @test NTL.roots(p)==[-1]
    @test NTL.integrate(Polynomials.derivative(p))==p-1
    for p in (2,3,5,7,11,13,17,19), n in 1:12
        c=zeros(UInt,n+1)
        ok=ccall((:_nmod_poly_conway,Nemo.libflint),Cint,(Ptr{UInt},UInt,Int),c,p,n)
        @test ok==1 && rep.(coefficients(conway_polynomial(p,n)))==c
    end
end
@testset "independent integer factor oracle" begin
    R,t=Nemo.polynomial_ring(Nemo.ZZ,:t)
    x=xpoly(BigInt)
    for a in -2:2, b in 1:3
        f=(a*x+1)*(x+b)^2*(x^2+1)
        z=NTL._factor_integer_kronecker(f)
        @test reconstruct(z)==f
        @test all(Nemo.is_irreducible(R(Nemo.ZZ.(coefficients(h)))) for (h,e) in z.factors)
        @test Dict(Tuple(coefficients(h))=>e for (h,e) in z.factors)==Dict(Tuple(coefficients(h))=>e for (h,e) in factor_integer_poly(f).factors)
    end
    @test ModPoly(Polynomials.Polynomial([1//2,1//3]))==poly(Rational{BigInt},[1//2,1//3])
end
@testset "integer kernel saturation and column lattice" begin
    for m in 1:4, n in 2:5, t in 1:3
        A=rand(rng,-7:7,m,n)
        t==1 && (A[:,1]=2A[:,2])
        K=kernel_integer(A)
        I=image_integer(A)
        h=Nemo.hnf(Nemo.matrix(Nemo.ZZ,permutedims(A)))
        canon=Nemo.hnf(Nemo.matrix(Nemo.ZZ,permutedims(I)))
        @test all(h[i,j]==canon[i,j] for i in axes(I,2),j in 1:m)
        @test all(iszero(h[i,j]) for i in size(I,2)+1:n,j in 1:m)
        if size(K,2)>0
            d=Nemo.snf(Nemo.matrix(Nemo.ZZ,K))
            @test all(isone(d[i,i]) for i in axes(K,2))
        end
    end
end
