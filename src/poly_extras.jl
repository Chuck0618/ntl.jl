export resultant, discriminant, subresultants, integrate, roots
Base.divrem(f::ModPoly,g::ModPoly)=DivRem(f,g)
roots(f::ModPoly)=FindRoots(f)
function integrate(f::ModPoly)
    F=field(f)
    result=[zero(F)]
    for (k,c) in enumerate(f._coeffs)
        if F===BigInt
            q,r=divrem(c,k)
            iszero(r) || throw(DomainError(c,"no integral coefficient antiderivative"))
            push!(result,q)
        elseif iszero(F(k))
            iszero(c) || throw(DomainError(k,"no antiderivative in this characteristic"))
            push!(result,zero(F))
        else
            push!(result,c/F(k))
        end
    end
    poly(F,result)
end
function _ringdet(M,F)
    n=size(M,1)
    c=ConstTerm(charpoly(M;ring=F))
    isodd(n) ? -c : c
end
"""Determinantal subresultants in increasing index, for indices 0 through min(m,n)-1."""
function subresultants(f::ModPoly,g::ModPoly)
    F=_check_poly_parents(f,g)
    m,n=degree(f),degree(g)
    out=typeof(f)[]
    for j in 0:min(m,n)-1
        h=m+n-2j
        S=fill(zero(F),h,h+j)
        for (row,(p,shift)) in enumerate(Iterators.flatten((((f,k) for k in n-j-1:-1:0),((g,k) for k in m-j-1:-1:0))))
            for d in 0:degree(p)
                S[row,h+j-d-shift]=coeff(p,d)
            end
        end
        push!(out,poly(F,[_ringdet(S[:,vcat(1:h-1,h+j-k)],F) for k in 0:j]))
    end
    out
end
function resultant(f::ModPoly,g::ModPoly)
    F=_check_poly_parents(f,g)
    (iszero(f)||iszero(g)) && return zero(F)
    degree(g)==0 && return ConstTerm(g)^degree(f)
    degree(f)==0 && return ConstTerm(f)^degree(g)
    ConstTerm(first(subresultants(f,g)))
end
function discriminant(f::ModPoly)
    n=degree(f)
    n>=0 || throw(DomainError(f))
    n<=1 && return one(field(f))
    r=resultant(f,derivative(f))
    # A characteristic-dependent degree drop of the derivative removes powers
    # of the leading coefficient from the Sylvester determinant.
    k=n-1-degree(derivative(f))
    r *= LeadCoeff(f)^k
    q=field(f)===BigInt ? div(r,LeadCoeff(f)) : r/LeadCoeff(f)
    isodd(n*(n-1)÷2) ? -q : q
end
