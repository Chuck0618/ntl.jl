export factor_integer_poly, content_primitive, is_irreducible_integer_poly
function content_primitive(f::ModPoly{BigInt})
    g = foldl(gcd,f._coeffs;init=big(0))
    iszero(g) && return (;unit=big(0),content=g,primitive=zero(f))
    u = sign(LeadCoeff(f))
    (;unit=u,content=g,primitive=poly(BigInt,[div(c,u*g) for c in f._coeffs]))
end

# Exact Kronecker interpolation; finite reference algorithm for small inputs.
function _signed_divisors(n)
    ds = BigInt[1]
    for (p,e) in FactorInt(abs(n))
        ds = [d*p^i for d in ds for i in 0:e]
    end
    vcat(ds,-ds)
end

function _integer_quotient(f,g)
    r = coefficients(f)
    q = zeros(BigInt,degree(f)-degree(g)+1)
    for k in reverse(eachindex(q))
        q[k],t = divrem(r[k+degree(g)],LeadCoeff(g))
        iszero(t) || return nothing
        for j in eachindex(g._coeffs)
            r[k+j-1] -= q[k]*g._coeffs[j]
        end
    end
    all(iszero,r) ? poly(BigInt,q) : nothing
end

function _split_integer(f)
    R = Rational{BigInt}
    for d in 1:degree(f)÷2
        xs = BigInt[]
        choices = Vector{Vector{BigInt}}()
        for j in 0:d
            x = big(j)
            y = evaluate(f,x)
            if iszero(y)
                g = poly(BigInt,[-x,1])
                return g,_integer_quotient(f,g)
            end
            push!(xs,x)
            push!(choices,_signed_divisors(y))
        end
        for ys in Iterators.product(choices...)
            h = interpolate(R,xs,collect(ys))
            degree(h)==d && all(a->denominator(a)==1,h._coeffs) || continue
            g = poly(BigInt,[numerator(a) for a in h._coeffs])
            LeadCoeff(g)>0 || continue
            q = _integer_quotient(f,g)
            q === nothing || return g,q
        end
    end
    nothing
end

function factor_integer_poly(f::ModPoly{BigInt})
    iszero(f) && throw(DomainError(f,"zero has no finite factorization"))
    ext = Base.get_extension(@__MODULE__, :NTLNemoExt)
    isnothing(ext) || return ext.factor_zz(f)
    _factor_integer_kronecker(f)
end

function _factor_integer_kronecker(f)
    u,c,g = content_primitive(f)
    pending = [g]
    factors = Tuple{typeof(f),Int}[]
    while !isempty(pending)
        h = pop!(pending)
        degree(h)==0 && continue
        parts = _split_integer(h)
        if parts === nothing
            idx = findfirst(t->first(t)==h, factors)
            if isnothing(idx)
                push!(factors,(h,1))
            else
                factors[idx] = (h,factors[idx][2]+1)
            end
        else
            append!(pending,parts)
        end
    end
    PolyFactorization(u*c,factors)
end
function is_irreducible_integer_poly(f::ModPoly{BigInt})
    degree(f)>0 || return false
    z = factor_integer_poly(f)
    abs(z.unit)==1 && length(z.factors)==1 && only(z.factors)[2]==1
end
function reconstruct(z::PolyFactorization{BigInt})
    foldl((p,t)->p*t[1]^t[2],z.factors;init=poly(BigInt,z.unit))
end
factor(f::ModPoly{BigInt})=factor_integer_poly(f)
is_irreducible(f::ModPoly{BigInt})=is_irreducible_integer_poly(f)
function gcd(f::ModPoly{BigInt},g::ModPoly{BigInt})
    Q=Rational{BigInt}
    h=gcd(poly(Q,coefficients(f)),poly(Q,coefficients(g)))
    iszero(h) && return zero(f)
    scale=foldl(lcm,denominator.(h._coeffs);init=big(1))
    c=gcd(content_primitive(f).content,content_primitive(g).content)
    c*content_primitive(poly(BigInt,[BigInt(x*scale) for x in h._coeffs])).primitive
end
function DivRem(f::ModPoly{BigInt},g::ModPoly{BigInt})
    Q=Rational{BigInt}
    q,r=DivRem(poly(Q,coefficients(f)),poly(Q,coefficients(g)))
    all(isinteger,q._coeffs) && all(isinteger,r._coeffs) || throw(DomainError(g,"division requires rational coefficients"))
    poly(BigInt,coefficients(q)),poly(BigInt,coefficients(r))
end
div(f::ModPoly{BigInt},g::ModPoly{BigInt})=first(DivRem(f,g))
rem(f::ModPoly{BigInt},g::ModPoly{BigInt})=last(DivRem(f,g))
function roots(f::ModPoly{BigInt})
    z=factor_integer_poly(f)
    out=BigInt[]
    for (h,e) in z.factors
        if degree(h)==1
            x,r=divrem(-ConstTerm(h),LeadCoeff(h))
            iszero(r) && append!(out,fill(x,e))
        end
    end
    out
end
