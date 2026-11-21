# Arbitrary-precision modular integers ----------------------------------------

export ZZPField, ZZ_p, ZZP, parent

"""A reusable parent object for arithmetic modulo an arbitrary positive BigInt."""
struct ZZPField
    _modulus::BigInt
    function ZZPField(n::Integer)
        n > 0 || throw(DomainError(n, "the modulus must be positive"))
        new(BigInt(n))
    end
end

"""Construct a reusable arbitrary-precision modular parent."""
ZZP(n::Integer) = ZZPField(n)

"""An arbitrary-precision modular integer. Use `F = ZZP(n); F(a)` to construct."""
struct ZZ_p <: Number
    _parent::ZZPField
    _rep::BigInt
    function ZZ_p(F::ZZPField, a::Integer)
        new(F, mod(BigInt(a), modulus(F)))
    end
end

(F::ZZPField)(a::Integer) = ZZ_p(F, a)
(F::ZZPField)(a::ZZ_p) = (parent(a) == F || throw(ArgumentError("modulus mismatch")); a)

modulus(F::ZZPField) = F._modulus
modulus(a::ZZ_p) = modulus(parent(a))
characteristic(F::ZZPField) = modulus(F)
characteristic(a::ZZ_p) = modulus(a)
parent(a::ZZ_p) = a._parent
parent(::zzp{P}) where {P} = zzp{P}
rep(a::ZZ_p) = a._rep

==(F::ZZPField, G::ZZPField) = modulus(F) == modulus(G)
isequal(F::ZZPField, G::ZZPField) = F == G
hash(F::ZZPField, h::UInt) = hash(modulus(F), h)

zero(F::ZZPField) = F(0)
one(F::ZZPField) = F(1)
oneunit(F::ZZPField) = F(1)
zero(a::ZZ_p) = zero(parent(a))
one(a::ZZ_p) = one(parent(a))
oneunit(a::ZZ_p) = one(a)
iszero(a::ZZ_p) = iszero(rep(a))
isone(a::ZZ_p) = isone(rep(a))

function show(io::IO, F::ZZPField)
    print(io, "ZZP(", modulus(F), ")")
end
function show(io::IO, a::ZZ_p)
    print(io, "ZZP(", modulus(a), ")(", rep(a), ")")
end

function _check_same_parent(x::ZZ_p, y::ZZ_p)
    parent(x) == parent(y) || throw(ArgumentError("cannot mix different moduli"))
    parent(x)
end

+(x::ZZ_p, y::ZZ_p) = (F = _check_same_parent(x, y); F(rep(x) + rep(y)))
-(x::ZZ_p, y::ZZ_p) = (F = _check_same_parent(x, y); F(rep(x) - rep(y)))
-(x::ZZ_p) = parent(x)(-rep(x))
*(x::ZZ_p, y::ZZ_p) = (F = _check_same_parent(x, y); F(rep(x) * rep(y)))

+(x::ZZ_p, y::Integer) = x + parent(x)(y)
+(x::Integer, y::ZZ_p) = parent(y)(x) + y
-(x::ZZ_p, y::Integer) = x - parent(x)(y)
-(x::Integer, y::ZZ_p) = parent(y)(x) - y
*(x::ZZ_p, y::Integer) = x * parent(x)(y)
*(x::Integer, y::ZZ_p) = parent(y)(x) * y

function inv(x::ZZ_p)
    n = modulus(x)
    gcd(rep(x), n) == 1 || throw(DomainError(x, "inverse undefined for a nonunit"))
    parent(x)(invmod(rep(x), n))
end
/(x::ZZ_p, y::ZZ_p) = x * inv(y)
/(x::ZZ_p, y::Integer) = x / parent(x)(y)
/(x::Integer, y::ZZ_p) = parent(y)(x) / y
sqr(x::ZZ_p) = x * x

function power(x::ZZ_p, e::Integer)
    base = x
    k = BigInt(e)
    if k < 0
        base = inv(base)
        k = -k
    end
    result = one(parent(x))
    while k > 0
        isodd(k) && (result *= base)
        k >>= 1
        k > 0 && (base *= base)
    end
    result
end
^(x::ZZ_p, e::Integer) = power(x, e)

==(x::ZZ_p, y::ZZ_p) = parent(x) == parent(y) && rep(x) == rep(y)
isequal(x::ZZ_p, y::ZZ_p) = x == y
==(x::ZZ_p, y::Integer) = rep(x) == mod(BigInt(y), modulus(x))
==(x::Integer, y::ZZ_p) = y == x
isequal(x::ZZ_p, y::Integer) = x == y
isequal(x::Integer, y::ZZ_p) = y == x
hash(x::ZZ_p, h::UInt) = hash(rep(x), h)
convert(::Type{BigInt}, x::ZZ_p) = BigInt(rep(x))
Base.BigInt(x::ZZ_p) = rep(x)

function _random_big_below(rng::AbstractRNG, n::BigInt)
    n > 0 || throw(DomainError(n, "the random bound must be positive"))
    bits = ndigits(n, base = 2)
    nbytes = cld(bits, 8)
    excess = 8 * nbytes - bits
    while true
        x = BigInt(0)
        for _ in 1:nbytes
            x = (x << 8) | BigInt(rand(rng, UInt8))
        end
        x >>= excess
        x < n && return x
    end
end

Random.rand(F::ZZPField) = Random.rand(Random.default_rng(), F)
Random.rand(rng::AbstractRNG, F::ZZPField) = F(_random_big_below(rng, modulus(F)))
Random.rand(F::ZZPField, dims::Dims) = Random.rand(Random.default_rng(), F, dims)
Random.rand(F::ZZPField, d::Integer, dims::Integer...) = Random.rand(Random.default_rng(), F, (d, dims...))
function Random.rand(rng::AbstractRNG, F::ZZPField, dims::Dims)
    A = Array{ZZ_p}(undef, dims)
    for i in eachindex(A)
        A[i] = rand(rng, F)
    end
    A
end
function Random.rand(rng::AbstractRNG, F::ZZPField, d::Integer, dims::Integer...)
    rand(rng, F, (d, dims...))
end
random_zzp(F::ZZPField; rng::AbstractRNG = Random.default_rng()) = rand(rng, F)
VectorRandom(F::ZZPField, k::Integer; rng::AbstractRNG = Random.default_rng()) = rand(rng, F, Int(k))
VectorConv(F::ZZPField, a::AbstractVector{<:Integer}) = [F(x) for x in a]
