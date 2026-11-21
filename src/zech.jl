# Zech/discrete-log representation for small finite extension fields ------------

export ZechField, ZechElem, zech_field, extension_element

struct ZechField{T,F}
    _extension::ExtensionField{T,F}
    _exp::Vector{Int}
    _log::Vector{Int}
    _generator::Int
    _powers::Vector{Int}
    _id::Int
    _p::Int
    _d::Int
end

struct ZechElem{F} <: Number
    _rep::Int
    _parent_id::Int
end

const _ZECH_REGISTRY = Dict{Int,Any}()
const _ZECH_NEXT_ID = Ref(0)
function _next_zech_id()
    _ZECH_NEXT_ID[] += 1
end

ZechElem(F::ZechField, v::Int) = ZechElem{typeof(F)}(v, F._id)
@inline parent(a::ZechElem{F}) where {F} = _ZECH_REGISTRY[a._parent_id]::F
@inline rep(a::ZechElem) = a._rep
modulus(F::ZechField) = modulus(F._extension)
modulus(a::ZechElem) = modulus(parent(a))
characteristic(F::ZechField) = characteristic(F._extension)
characteristic(a::ZechElem) = characteristic(parent(a))
cardinality(F::ZechField) = length(F._exp) + 1
extension_degree(F::ZechField) = extension_degree(F._extension)
extension_field(F::ZechField) = F._extension

function _v_of(a::ExtensionElem)
    f = rep(a)
    p = Int(characteristic(a))
    v = 0
    for i in 0:degree(f)
        v += Int(rep(coeff(f, i))) * p^i
    end
    v
end

function _poly_from_v(F::Type{zzp{P}}, v::Int, d::Int) where {P}
    c = Int[]
    for _ in 1:d
        push!(c, mod(v, P))
        v ÷= P
    end
    poly(F, c)
end
_poly_from_v(F, v::Int, d::Int) = poly(F, [mod(v, modulus(F)) for _ in 1:d])

function _extension_elements(E::ExtensionField)
    F = base_field(E)
    d = extension_degree(E)
    (E(collect(c)) for c in Iterators.product(ntuple(_ -> F, d)...))
end

function _multiplicative_order(a::ExtensionElem, q::Int)
    y = one(a)
    for k in 1:(q - 1)
        y *= a
        y == one(a) && return k
    end
    q
end

function _find_generator(E::ExtensionField, q::Int)
    factors = FactorInt(q - 1)
    for a in _extension_elements(E)
        iszero(a) && continue
        isone(a) && continue
        if all(a^div(q - 1, r) != one(a) for (r, _) in factors)
            return a
        end
    end
    error("no primitive element found")
end

function _check_generator(a::ExtensionElem, q::Int)
    iszero(a) && throw(ArgumentError("zero cannot be a multiplicative generator"))
    _multiplicative_order(a, q) == q - 1 || throw(ArgumentError("the supplied element is not primitive"))
    factors = FactorInt(q - 1)
    all(a^div(q - 1, r) != one(a) for (r, _) in factors) || throw(ArgumentError("the supplied element is not primitive"))
    a
end

"""Build a Zech representation for an extension of order at most 65535."""
function zech_field(E::ExtensionField; primitive = nothing)
    F = base_field(E)
    F isa Type && F <: zzp || throw(ArgumentError("Zech representation currently requires a machine-word prime base field"))
    q = Int(cardinality(E))
    q <= 65535 || throw(ArgumentError("Zech tables are limited to fields with at most 65535 elements"))
    g = primitive === nothing ? _find_generator(E, q) : _check_generator(primitive, q)
    gv = _v_of(g)
    exp = Int[]
    log = zeros(Int, q)
    y = one(E)
    for k in 0:(q - 2)
        v = _v_of(y)
        push!(exp, v)
        log[v + 1] = k
        y *= g
    end
    y == one(E) || error("internal error while building the Zech table")
    id = _next_zech_id()
    F = ZechField(E, exp, log, gv, [Int(characteristic(E))^i for i in 0:(extension_degree(E) - 1)], id, Int(characteristic(E)), extension_degree(E))
    _ZECH_REGISTRY[id] = F
    F
end

@inline zero(F::ZechField) = ZechElem(F, 0)
@inline one(F::ZechField) = ZechElem(F, 1)
@inline zero(a::ZechElem) = zero(parent(a))
@inline one(a::ZechElem) = one(parent(a))
oneunit(a::ZechElem) = one(a)
iszero(a::ZechElem) = iszero(rep(a))
isone(a::ZechElem) = isone(rep(a))
==(a::ZechElem, b::ZechElem) = parent(a) == parent(b) && rep(a) == rep(b)
isequal(a::ZechElem, b::ZechElem) = a == b
hash(a::ZechElem, h::UInt) = hash(rep(a), h)
Base.Int(a::ZechElem) = rep(a)
Base.length(F::ZechField) = cardinality(F)
function Base.iterate(F::ZechField, state::Int = 0)
    state >= cardinality(F) && return nothing
    (ZechElem(F, state), state + 1)
end

"""Return the ordinary polynomial extension element represented by a Zech element."""
extension_element(a::ZechElem) = parent(a)._extension(_poly_from_v(base_field(parent(a)._extension), rep(a), extension_degree(parent(a))))

function show(io::IO, a::ZechElem)
    show(io, extension_element(a))
end

@inline function _check_zech_parents(a::ZechElem, b::ZechElem)
    a._parent_id === b._parent_id || throw(ArgumentError("different extension parents"))
    parent(a)
end

@inline function _add_v(F::ZechField, u::Int, v::Int)
    p = F._p
    d = F._d
    if d == 2
        return mod(u + v, p) + p * mod((u ÷ p) + (v ÷ p), p)
    end
    w = 0
    @inbounds for pow in F._powers
        w += mod((u ÷ pow) + (v ÷ pow), p) * pow
    end
    w
end

@inline function _sub_v(F::ZechField, u::Int, v::Int)
    p = F._p
    d = F._d
    if d == 2
        return mod(u - v, p) + p * mod((u ÷ p) - (v ÷ p), p)
    end
    w = 0
    @inbounds for pow in F._powers
        w += mod((u ÷ pow) - (v ÷ pow), p) * pow
    end
    w
end

@inline function _mul_v(F::ZechField, u::Int, v::Int)
    (iszero(u) || iszero(v)) && return 0
    F._exp[mod(F._log[u + 1] + F._log[v + 1], length(F._exp)) + 1]
end

@inline function _inv_v(F::ZechField, u::Int)
    iszero(u) && throw(DomainError(u, "zero has no inverse"))
    F._exp[mod(-F._log[u + 1], length(F._exp)) + 1]
end

@inline function _pow_v(F::ZechField, u::Int, e::Int)
    if iszero(u)
        e == 0 && return 1
        e < 0 && throw(DomainError(u, "zero cannot be raised to a negative power"))
        return 0
    end
    F._exp[mod(Int128(e) * F._log[u + 1], length(F._exp)) + 1]
end

@inline function _pow_v(F::ZechField, u::Int, e::Integer)
    if iszero(u)
        e == 0 && return 1
        e < 0 && throw(DomainError(u, "zero cannot be raised to a negative power"))
        return 0
    end
    F._exp[mod(BigInt(e) * F._log[u + 1], length(F._exp)) + 1]
end

@inline +(a::ZechElem, b::ZechElem) = (F = _check_zech_parents(a, b); ZechElem(F, _add_v(F, rep(a), rep(b))))
@inline -(a::ZechElem, b::ZechElem) = (F = _check_zech_parents(a, b); ZechElem(F, _sub_v(F, rep(a), rep(b))))
-(a::ZechElem) = (F = parent(a); iszero(a) ? a : ZechElem(F, _sub_v(F, 0, rep(a))))
@inline *(a::ZechElem, b::ZechElem) = (F = _check_zech_parents(a, b); ZechElem(F, _mul_v(F, rep(a), rep(b))))
@inline /(a::ZechElem, b::ZechElem) = a * inv(b)
@inline @inline muladd(x::ZechElem, y::ZechElem, z::ZechElem) = begin
    F = _check_zech_parents(x, y)
    _check_zech_parents(x, z)
    v = _mul_v(F, rep(x), rep(y))
    ZechElem(F, _add_v(F, v, rep(z)))
end

@inline inv(a::ZechElem) = ZechElem(parent(a), _inv_v(parent(a), rep(a)))
@inline ^(a::ZechElem, e::Integer) = ZechElem(parent(a), _pow_v(parent(a), rep(a), e))
sqr(a::ZechElem) = a * a
power(a::ZechElem, e::Integer) = a^e

Random.rand(F::ZechField) = ZechElem(F, rand(0:(cardinality(F) - 1)))
Random.rand(rng::AbstractRNG, F::ZechField) = ZechElem(F, rand(rng, 0:(cardinality(F) - 1)))
function Random.rand(rng::AbstractRNG, F::ZechField, dims::Dims)
    A = Array{ZechElem}(undef, dims)
    for i in eachindex(A)
        A[i] = rand(rng, F)
    end
    A
end
Random.rand(F::ZechField, dims::Dims) = rand(Random.default_rng(), F, dims)
Random.rand(rng::AbstractRNG, F::ZechField, d::Integer, dims::Integer...) = rand(rng, F, (d, dims...))
Random.rand(F::ZechField, d::Integer, dims::Integer...) = rand(Random.default_rng(), F, (d, dims...))
