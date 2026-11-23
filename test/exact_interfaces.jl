using Test, NTL, Random
@testset "matrix and quotient polynomials" begin
    rng=MersenneTwister(444)
    E,a=conway_field(2,3)
    for F in (BigInt,GF2,PrimeField(5),BigModRing(12),E)
        for n in 1:4
            M=[F(rand(rng,-4:4)) for _ in 1:n,_ in 1:n]
            F===E && (M[1,1]=a)
            f=charpoly(M)
            @test degree(f)==n && isone(LeadCoeff(f))
            e=fill(zero(F),n,n)
            for c in reverse(coefficients(f))
                e=[sum((M[i,k]*e[k,j] for k in 1:n);init=zero(F))+(i==j ? c : zero(F)) for i in 1:n,j in 1:n]
            end
            @test all(iszero,e)
        end
        @test charpoly(fill(zero(F),3,3))==xpoly(F)^3
        @test charpoly([i==j ? one(F) : zero(F) for i=1:3,j=1:3])==(xpoly(F)-one(F))^3
    end
    @test_throws DimensionMismatch charpoly(zeros(Int,2,3))
    @test charpoly(zeros(Int,0,0))==poly(BigInt,1)
    for bits in 0:7, hbits in 0:7
        f=to_modpoly(GF2X(8+bits))
        h=to_modpoly(GF2X(hbits))
        g=MinPolyMod(h,f)
        @test iszero(CompMod(g,h,f)) && ismonic(g)
        for i in 1:(1<<degree(g))-1
            @test !iszero(CompMod(to_modpoly(GF2X(i)),h,f))
        end
    end
    @test minpoly(a)==defining_polynomial(E)==charpoly(a)
    @test charpoly(one(E))==(xpoly(GF2)-1)^3
    @test_throws ArgumentError MinPolyMod(xpoly(BigModRing(6)),xpoly(BigModRing(6))^2)
end
@testset "integer lattices" begin
    rng=MersenneTwister(104)
    for m in 0:5,n in 0:5,t in 1:4
        A=rand(rng,-5:5,m,n)
        t==2 && fill!(A,0)
        H,U=hnf(A)
        @test U*A==H
        m>0 && (@test abs(det_exact(U))==1)
        @test hnf(H;transform=false)==H
        piv=[findfirst(!iszero,H[i,:]) for i in 1:m]
        non=filter(!isnothing,piv)
        @test issorted(non;lt=(<)) && allunique(non)
        @test all(isnothing,piv[length(non)+1:end])
        for (i,j) in enumerate(non)
            @test H[i,j]>0 && all(x->0<=x<H[i,j],H[1:i-1,j])
        end
        K=kernel_integer(A)
        @test iszero(A*K)
        @test size(image_integer(A),2)==length(non)==n-size(K,2)
        b=A*rand(rng,-3:3,n)
        y=solve_integer(A,b)
        @test y!==nothing && A*y.x==b && K==y.kernel
    end
    @test solve_integer([2 4],[1])===nothing
    @test image_integer([6 10])==reshape(BigInt[2],1,1)
    @test hnf([big(2)^140 1;1 0]).H==[1 0;0 1]
end
@testset "factorization and domains" begin
    x=xpoly(BigInt)
    for f in (-6*(2*x+1)^2*x*(x^2+1),poly(BigInt,-9),x^4+1,(x+big(2)^80)*x)
        z=factor_integer_poly(f)
        @test reconstruct(z)==f
        @test all(t->is_irreducible_integer_poly(t[1]),z.factors)
    end
    @test is_irreducible_integer_poly(x^4+1)
    @test !is_irreducible_integer_poly(x^4-1)
    @test roots(2*x^2-x-1)==[1]
    @test content_primitive(-6*x).primitive==x
    @test_throws DomainError factor_integer_poly(zero(x))
    for F in (GF2,PrimeField(3),BigModRing(big(101)))
        E=extension_field(F,BuildIrred(F,2;rng=MersenneTwister(55)))
        a=gen(E); x=xpoly(E)
        f=(x+a)^2*(x+a+1)
        z=characteristic(E)==2 ? GF2EXFactoring(f) : ZZ_pEXFactoring(f)
        @test reconstruct(z)==f && (x+a,2) in z.factors
        @test minpoly(a)==defining_polynomial(E)
    end
    E,a=conway_field(2,3;symbol=:beta)
    @test occursin("beta", sprint(show,E))
    T=GF2E(GF2X(13))
    phi=identify(E,T)
    psi=identify(T,E;generator_map=first(b for b in NTL._field_values(E) if phi(b)==gen(T)))
    for v in NTL._field_values(E)
        @test psi(phi(v))==v
        @test phi(v+a)==phi(v)+phi(a) && phi(v*a)==phi(v)*phi(a)
        @test trace(v)==trace(phi(v)) && norm(v)==norm(phi(v))
        !iszero(v) && (@test exp(log(v),E)==v && isone(v^order(v)))
    end
    @test_throws ArgumentError identify(E,T;generator_map=one(T))
    @test_throws KeyError conway_polynomial(2,10000)
    Z=enable_zech_multiplication(E)
    for u in Z,v in Z
        @test extension_element(u*v)==extension_element(u)*extension_element(v)
        @test Z(extension_element(u))==u
    end
    @test disable_zech_multiplication(Z)==E
    @test_throws DomainError log(zero(Z))
end
@testset "storage, polynomial extras, BKZ, real arithmetic and tasks" begin
    z=ZZVec([-big(2)^75,2,0])
    @test collect(z+z)==2 .* collect(z)
    @test copy(z) !== z && collect(z[1:2])==[-big(2)^75,2]
    @test packed_dot(z,z)==big(2)^150+4
    @test isempty(ZZVec(Int[]))
    @test_throws BoundsError z[0]
    @test_throws InexactError setindex!(z,big(2)^999,1)
    g=GF2XVec(GF2X.([3,6]))
    @test packed_dot(g,g)==GF2X(3)^2+GF2X(6)^2
    for (p,T) in ((97,Int8),(257,Int16),(65537,Int32),(big(2)^63+1,Int128),(big(2)^180,BigInt))
        R=BigModRing(p)
        v=packedfieldvec(R,[0,p-1,-1])
        @test eltype(v.data)==T && v[2]==R(-1) && v[3]==R(p-1)
    end
    x=xpoly(BigInt)
    @test resultant(x^2+1,x+1)==2
    @test discriminant(3*x^2+4*x-1)==28
    @test derivative(integrate(6*x^2-2*x+1))==6*x^2-2*x+1
    @test length(isolate_real_roots((x-1)^2*(x^2-2)))==3
    for A in ([1 1 0;2 0 1;3 1 2],[2 3 1 4;0 1 1 2;1 2 1 4])
        z=bkz(A;block_size=3)
        @test z.verified && z.transform*A==z.basis && abs(det_exact(z.transform))==1
        # Independent bounded brute force on the first (unprojected) full block.
        bound=sum(abs2,z.basis[1,:])
        for c in Iterators.product(ntuple(_->-3:3,3)...)
            all(iszero,c) && continue
            @test sum(abs2,permutedims(collect(c))*z.basis)>=bound
        end
    end
    @test abs(Float64(sin(quad_float(1)))-sin(1))<1e-15
    @test sqrt(xdouble(4))==xdouble(2)
    table=LazyTable{Int,Int}()
    count=Threads.Atomic{Int}(0)
    values=parallel_map(BasicThreadPool(4),1:20) do i
        get!(table,1) do
            Threads.atomic_add!(count,1)
            19
        end
    end
    @test values==fill(19,20) && count[]==1
    @test Lazy(()->42)[]==42
end
@testset "parallel caches and arithmetic boundaries" begin
    c=ArithmeticCache()
    F=cached_field(c,97)
    a=F.([1,2,3]); b=F.([5,6])
    @test all(v->v==convolve(a,b),parallel_map(_->cached_ntt_convolve(c,a,b),BasicThreadPool(4),1:20))
    E,g=conway_field(2,3)
    @test cached_zech(c,E)===cached_zech(c,E)
    @test length(c.roots.values)==1 && length(c.twiddles.values)==2
    @test order(one(GF2))==1 && is_primitive_root(one(GF2))
    @test_throws ArgumentError bkz([1 1;2 2])
    @test order(PrimeField(7)(2))==3
    @test Base.exp(log(PrimeField(7)(3)),PrimeField(7))==PrimeField(7)(3)
    @test_throws ArgumentError identify(E,conway_field(2,2)[1])
end
@testset "BKZ insertion and independent projected-block certificate" begin
    A=[28 18 -1 -6;29 21 -7 27;-21 6 11 5;-28 -9 24 -7]
    B=lll(A;delta=99//100,mode=:exact).basis
    result=bkz(A;block_size=4)
    @test sum(abs2,result.basis[1,:])<sum(abs2,B[1,:])
    @test result.transform*A==result.basis && abs(det_exact(result.transform))==1
    G=Rational{BigInt}.(result.basis*permutedims(result.basis))
    for lo in 1:3
        I=lo:4
        K=G[I,I]
        if lo>1
            J=1:lo-1
            K-=G[I,J]*matrix_inverse(G[J,J])*G[J,I]
        end
        D=matrix_inverse(K)
        # Dual Cauchy–Schwarz bounds every coordinate of a shorter vector.
        limits=[isqrt(floor(BigInt,K[1,1]*D[i,i])) for i in axes(K,1)]
        for v in Iterators.product(((-a:a) for a in limits)...)
            all(iszero,v) && continue
            z=collect(v)
            @test sum(z.*(K*z))>=K[1,1]
        end
    end
end
