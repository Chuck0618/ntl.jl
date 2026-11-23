export conway_field, conway_polynomial, identify, primitive_element, order, is_primitive_root
export any_primitive_root, minimal_primitive_root, enable_zech_multiplication, disable_zech_multiplication
export GF2EXFactoring, ZZ_pEXFactoring

function conway_polynomial(p::Integer,n::Integer)
    path = joinpath(@__DIR__,"..","data","conway.txt")
    open(path) do io
        for line in eachline(io)
            row = parse.(Int,split(line))
            (row[1],row[2]) == (p,n) && return poly(p==2 ? GF2 : PrimeField(p),row[3:end])
        end
        throw(KeyError((p,n)))
    end
end
function conway_field(p,n;symbol::Symbol=:α)
    m = conway_polynomial(p,n)
    E = ExtensionField(field(m),m,symbol)
    E,gen(E)
end
function _field_values(E::ExtensionField)
    p = characteristic(E)
    cardinality(base_field(E))==p || throw(ArgumentError("expected a prime base"))
    (E([div(k,p^i)%p for i in 0:degree(E)-1]) for k in big(0):cardinality(E)-1)
end
function identify(S::ExtensionField,T::ExtensionField;generator_map=nothing)
    characteristic(S)==characteristic(T) && degree(S)==degree(T) || throw(ArgumentError("incompatible extensions"))
    cardinality(base_field(S))==characteristic(S) && cardinality(base_field(T))==characteristic(T) || throw(ArgumentError("expected prime bases"))
    f = defining_polynomial(S)
    eval_at(h,a) = foldr((b,v)->v*a+T(rep(b)),h._coeffs;init=zero(T))
    root = generator_map
    if isnothing(root)
        if S==T
            root = gen(T)
        else
            for a in _field_values(T)
                if iszero(eval_at(f,a))
                    root = a
                    break
                end
            end
        end
    end
    root isa ExtensionElem && parent(root)==T && iszero(eval_at(f,root)) || throw(ArgumentError("invalid generator image"))
    function transport(a::ExtensionElem)
        parent(a)==S || throw(ArgumentError("wrong source"))
        eval_at(rep(a),root)
    end
end
Base.convert(E::ExtensionField,a::ExtensionElem;generator_map=nothing) = identify(parent(a),E;generator_map)(a)
const GroupElem = Union{ExtensionElem,ZechElem,PrimeField,GF2,ModInt}
function order(a::GroupElem)
    iszero(a) && throw(DomainError(a,"zero has no multiplicative order"))
    F = parent(a)
    _validate_factor_parent(F)
    n = big(cardinality(F))-1
    for (p,e) in FactorInt(n)
        for _ in 1:e
            isone(a^div(n,p)) || break
            n = div(n,p)
        end
    end
    n
end
is_primitive_root(a::GroupElem) = !iszero(a) && order(a)==cardinality(parent(a))-1
function minimal_primitive_root(F)
    _validate_factor_parent(F)
    values = F isa ExtensionField ? _field_values(F) : (F(i) for i in big(1):cardinality(F)-1)
    for a in values
        is_primitive_root(a) && return a
    end
    throw(ArgumentError("no primitive element"))
end
minimal_primitive_root(F::ZechField) = first(a for a in F if is_primitive_root(a))
primitive_element(F) = minimal_primitive_root(F)
primitive_element(F::ZechField) = ZechElem(F,F._generator)
any_primitive_root(F) = primitive_element(F)
Base.exp(k::Integer,F::Union{Type{<:PrimeField},Type{GF2},BigModRing,ExtensionField,ZechField}) = primitive_element(F)^mod(k,cardinality(F)-1)
function Base.log(a::GroupElem)
    iszero(a) && throw(DomainError(a,"zero has no discrete logarithm"))
    F = parent(a)
    g = primitive_element(F)
    v = one(a)
    for k in big(0):big(cardinality(F))-2
        v==a && return k
        v *= g
    end
    error("element outside generated group")
end
Base.log(a::ZechElem) = iszero(a) ? throw(DomainError(a)) : parent(a)._log[rep(a)+1]
enable_zech_multiplication(E::ExtensionField) = zech_field(E;primitive=primitive_element(E))
disable_zech_multiplication(Z::ZechField) = extension_field(Z)
function (Z::ZechField)(a::ExtensionElem)
    parent(a)==extension_field(Z) || throw(ArgumentError("incompatible parents"))
    ZechElem(Z,_v_of(a))
end
# Repeated linear factors over small extension fields are cheap to remove.
# The remaining positive-degree component still uses the generic algorithm.
function _extension_factor(f::ModPoly{<:ExtensionElem};rng=Random.default_rng())
    E=field(f)
    if cardinality(E)>32 || degree(f)>12 || degree(f)<=0 || cardinality(base_field(E))!=characteristic(E)
        return factor(f;rng)
    end
    u=LeadCoeff(f)
    g=monic(f)
    factors=Tuple{typeof(f),Int}[]
    for a in _field_values(E)
        degree(g)==0 && break
        h=poly(E,[-a,one(E)])
        count=0
        while degree(g)>0 && iszero(evaluate(g,a))
            g=div(g,h)
            count+=1
        end
        count>0 && push!(factors,(h,count))
    end
    if degree(g)>0
        append!(factors,factor(g;rng).factors)
    end
    PolyFactorization(u,factors)
end
function GF2EXFactoring(f::ModPoly{<:ExtensionElem};kw...)
    characteristic(field(f))==2 || throw(ArgumentError("expected characteristic two"))
    _extension_factor(f;kw...)
end
ZZ_pEXFactoring(f::ModPoly{<:ExtensionElem};kw...) = _extension_factor(f;kw...)
