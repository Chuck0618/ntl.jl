# Finite-field extensions ------------------------------------------------------

export ExtensionField, ExtensionElem, ZZpE, extension_field, GF2E
export base_field, defining_polynomial, extension_degree, cardinality, gen
export trace, norm

"""A finite-field parent `F[x]/(m)` with an irreducible defining polynomial."""
struct ExtensionField{T,F}
    _base::F
    _modulus::ModPoly{T,F}
end

# Internal marker: only the parent call may construct an already reduced rep.
struct _UnsafeExtensionRep end
const _unsafe_extension_rep = _UnsafeExtensionRep()

"""
    ZZpE{T,F} <: Number

An element of an `ExtensionField`, represented by a polynomial reduced modulo
the defining polynomial. Use `E = extension_field(F, m)` to build the parent,
then normally construct elements with `E([a0, a1, ...])`, `E(a)`, or `E(p)`.
`ExtensionElem` is the descriptive alias for this type; `ZZpE` is the
NTL-facing spelling. Coefficient vectors run from low degree to high degree.
"""
struct ZZpE{T,F} <: Number
    _parent::ExtensionField{T,F}
    _rep::ModPoly{T,F}
    function ZZpE{T,F}(parent::ExtensionField{T,F}, rep::ModPoly{T,F},
                       ::_UnsafeExtensionRep) where {T,F}
        new{T,F}(parent, rep)
    end
end

@inline _raw_extension(E::ExtensionField{T,F}, p::ModPoly{T,F}) where {T,F} =
    ZZpE{T,F}(E, p, _unsafe_extension_rep)

"""Build an extension parent and optionally verify that `m` is irreducible."""
function extension_field(F, m::ModPoly; check::Bool = true)
    field(m) == F || throw(ArgumentError("the defining polynomial has the wrong base field"))
    degree(m) > 0 || throw(DomainError(m, "the defining polynomial must have positive degree"))
    check && !is_irreducible(m) && throw(ArgumentError("the defining polynomial is not irreducible"))
    ExtensionField(F, monic(m))
end
extension_field(F, m::GF2X; check::Bool = true) = extension_field(F, to_modpoly(m); check = check)

"""Construct a binary extension field from a `GF2X` or a `ModPoly` over `GF2`."""
GF2E(m::Union{GF2X,ModPoly}; check::Bool = true) = extension_field(GF2, m; check = check)

base_field(E::ExtensionField) = E._base
defining_polynomial(E::ExtensionField) = E._modulus
extension_degree(E::ExtensionField) = degree(defining_polynomial(E))
degree(E::ExtensionField) = extension_degree(E)
modulus(E::ExtensionField) = defining_polynomial(E)
characteristic(E::ExtensionField) = characteristic(base_field(E))

cardinality(::Type{zzp{P}}) where {P} = BigInt(P)
cardinality(F::ZZPField) = modulus(F)
cardinality(::Type{GF2}) = BigInt(2)
cardinality(E::ExtensionField) = cardinality(base_field(E))^extension_degree(E)
_field_order(E::ExtensionField) = cardinality(E)
_field_characteristic(E::ExtensionField) = BigInt(characteristic(E))

==(E::ExtensionField, G::ExtensionField) = base_field(E) == base_field(G) && defining_polynomial(E) == defining_polynomial(G)
isequal(E::ExtensionField, G::ExtensionField) = E == G

function (E::ExtensionField)(p::ModPoly)
    field(p) == base_field(E) || throw(ArgumentError("coefficient parent mismatch"))
    _raw_extension(E, rem(p, defining_polynomial(E)))
end
(E::ExtensionField)(c::AbstractVector) = E(poly(base_field(E), c))
(E::ExtensionField)(c::Number) = E(poly(base_field(E), c))
(E::ExtensionField)(a::ZZpE) = parent(a) == E ? a : throw(ArgumentError("extension parent mismatch"))
ZZpE(E::ExtensionField, p::ModPoly) = E(p)
ZZpE(E::ExtensionField, c::AbstractVector) = E(c)
ZZpE(E::ExtensionField, c::Number) = E(c)
ZZpE(E::ExtensionField, a::ZZpE) = E(a)

"""Descriptive alias for `ZZpE`; ordinary construction can simply use `E(...)`."""
const ExtensionElem = ZZpE

parent(a::ZZpE) = a._parent
rep(a::ZZpE) = a._rep
base_field(a::ZZpE) = base_field(parent(a))
modulus(a::ZZpE) = defining_polynomial(parent(a))
characteristic(a::ZZpE) = characteristic(parent(a))
zero(E::ExtensionField) = E(poly(base_field(E), 0))
one(E::ExtensionField) = E(poly(base_field(E), 1))
oneunit(E::ExtensionField) = one(E)
zero(a::ZZpE) = zero(parent(a))
one(a::ZZpE) = one(parent(a))
oneunit(a::ZZpE) = one(a)
iszero(a::ZZpE) = iszero(rep(a))
isone(a::ZZpE) = isone(rep(a))

"""Return the image of the polynomial indeterminate in the extension."""
gen(E::ExtensionField) = E(xpoly(base_field(E)))

function show(io::IO, E::ExtensionField)
    print(io, "ExtensionField(", base_field(E), ", ", defining_polynomial(E), ")")
end
show(io::IO, a::ZZpE) = print(io, "ZZpE(", parent(a), ", ", rep(a), ")")

function _check_extension_parents(a::ZZpE, b::ZZpE)
    parent(a) == parent(b) || throw(ArgumentError("extension parents differ"))
    parent(a)
end
+(a::ZZpE, b::ZZpE) = (_check_extension_parents(a, b); parent(a)(rep(a) + rep(b)))
-(a::ZZpE, b::ZZpE) = (_check_extension_parents(a, b); parent(a)(rep(a) - rep(b)))
-(a::ZZpE) = parent(a)(-rep(a))
*(a::ZZpE, b::ZZpE) = (_check_extension_parents(a, b); parent(a)(MulMod(rep(a), rep(b), defining_polynomial(parent(a)))))
sqr(a::ZZpE) = a * a

+(a::ZZpE, c::Number) = a + parent(a)(c)
+(c::Number, a::ZZpE) = parent(a)(c) + a
-(a::ZZpE, c::Number) = a - parent(a)(c)
-(c::Number, a::ZZpE) = parent(a)(c) - a
*(a::ZZpE, c::Number) = a * parent(a)(c)
*(c::Number, a::ZZpE) = parent(a)(c) * a

function inv(a::ZZpE)
    iszero(a) && throw(DomainError(a, "zero is not invertible"))
    parent(a)(InvMod(rep(a), defining_polynomial(parent(a))))
end
/(a::ZZpE, b::ZZpE) = a * inv(b)
/(a::ZZpE, c::Number) = a / parent(a)(c)
/(c::Number, a::ZZpE) = parent(a)(c) / a

function power(a::ZZpE, e::Integer)
    parent(a)(PowerMod(rep(a), e, defining_polynomial(parent(a))))
end
^(a::ZZpE, e::Integer) = power(a, e)

==(a::ZZpE, b::ZZpE) = parent(a) == parent(b) && rep(a) == rep(b)
isequal(a::ZZpE, b::ZZpE) = a == b
==(a::ZZpE, c::Number) = a == parent(a)(c)
==(c::Number, a::ZZpE) = a == c

Random.rand(E::ExtensionField) = rand(Random.default_rng(), E)
function Random.rand(rng::AbstractRNG, E::ExtensionField)
    E(random_poly(base_field(E), extension_degree(E); rng = rng))
end
function Random.rand(rng::AbstractRNG, E::ExtensionField, dims::Dims)
    T = typeof(zero(E))
    A = Array{T}(undef, dims)
    for i in eachindex(A)
        A[i] = rand(rng, E)
    end
    A
end
Random.rand(rng::AbstractRNG, E::ExtensionField, d::Integer, dims::Integer...) = rand(rng, E, (d, dims...))
Random.rand(E::ExtensionField, dims::Dims) = rand(Random.default_rng(), E, dims)
Random.rand(E::ExtensionField, d::Integer, dims::Integer...) = rand(Random.default_rng(), E, (d, dims...))

"""Absolute trace from an extension to its base field."""
function trace(a::ZZpE)
    E = parent(a)
    q = cardinality(base_field(E))
    d = extension_degree(E)
    term = a
    result = zero(a)
    for _ in 1:d
        result += term
        term = term^q
    end
    degree(rep(result)) <= 0 || error("trace did not land in the base field")
    ConstTerm(rep(result))
end

"""Absolute norm from an extension to its base field."""
function norm(a::ZZpE)
    E = parent(a)
    q = cardinality(base_field(E))
    d = extension_degree(E)
    result = a^div(q^d - 1, q - 1)
    degree(rep(result)) <= 0 || error("norm did not land in the base field")
    ConstTerm(rep(result))
end
