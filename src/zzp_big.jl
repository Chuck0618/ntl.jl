# Arbitrary-precision modular integers ----------------------------------------

export ZZPField, ZZp, ZZ_p, ZZP, parent

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

"""An arbitrary-precision modular integer. Use `R = ZZP(n); R(a)` or `ZZp(R, a)`."""
struct ZZp <: Number
    _parent::ZZPField
    _rep::BigInt
    function ZZp(F::ZZPField, a::Integer)
        new(F, mod(BigInt(a), modulus(F)))
    end
end

const ZZ_p = ZZp

(F::ZZPField)(a::Integer) = ZZp(F, a)
(F::ZZPField)(a::ZZp) = (parent(a) == F || throw(ArgumentError("modulus mismatch")); a)

modulus(F::ZZPField) = F._modulus
modulus(a::ZZp) = modulus(parent(a))
characteristic(F::ZZPField) = modulus(F)
characteristic(a::ZZp) = modulus(a)
parent(a::ZZp) = a._parent
parent(::zzp{P}) where {P} = zzp{P}
rep(a::ZZp) = a._rep

==(F::ZZPField, G::ZZPField) = modulus(F) == modulus(G)
isequal(F::ZZPField, G::ZZPField) = F == G
hash(F::ZZPField, h::UInt) = hash(modulus(F), h)

zero(F::ZZPField) = F(0)
one(F::ZZPField) = F(1)
oneunit(F::ZZPField) = F(1)
zero(a::ZZp) = zero(parent(a))
one(a::ZZp) = one(parent(a))
oneunit(a::ZZp) = one(a)
iszero(a::ZZp) = iszero(rep(a))
isone(a::ZZp) = isone(rep(a))

function show(io::IO, F::ZZPField)
    print(io, "ZZp(", modulus(F), ")")
end
function show(io::IO, a::ZZp)
    print(io, "ZZp(", modulus(a), ")(", rep(a), ")")
end

function _check_same_parent(x::ZZp, y::ZZp)
    parent(x) == parent(y) || throw(ArgumentError("cannot mix different moduli"))
    parent(x)
end

+(x::ZZp, y::ZZp) = (F = _check_same_parent(x, y); F(rep(x) + rep(y)))
-(x::ZZp, y::ZZp) = (F = _check_same_parent(x, y); F(rep(x) - rep(y)))
-(x::ZZp) = parent(x)(-rep(x))
*(x::ZZp, y::ZZp) = (F = _check_same_parent(x, y); F(rep(x) * rep(y)))

+(x::ZZp, y::Integer) = x + parent(x)(y)
+(x::Integer, y::ZZp) = parent(y)(x) + y
-(x::ZZp, y::Integer) = x - parent(x)(y)
-(x::Integer, y::ZZp) = parent(y)(x) - y
*(x::ZZp, y::Integer) = x * parent(x)(y)
*(x::Integer, y::ZZp) = parent(y)(x) * y

function inv(x::ZZp)
    n = modulus(x)
    gcd(rep(x), n) == 1 || throw(DomainError(x, "inverse undefined for a nonunit"))
    parent(x)(invmod(rep(x), n))
end
/(x::ZZp, y::ZZp) = x * inv(y)
/(x::ZZp, y::Integer) = x / parent(x)(y)
/(x::Integer, y::ZZp) = parent(y)(x) / y
sqr(x::ZZp) = x * x

function power(x::ZZp, e::Integer)
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
^(x::ZZp, e::Integer) = power(x, e)

==(x::ZZp, y::ZZp) = parent(x) == parent(y) && rep(x) == rep(y)
isequal(x::ZZp, y::ZZp) = x == y
==(x::ZZp, y::Integer) = rep(x) == mod(BigInt(y), modulus(x))
==(x::Integer, y::ZZp) = y == x
isequal(x::ZZp, y::Integer) = x == y
isequal(x::Integer, y::ZZp) = y == x
hash(x::ZZp, h::UInt) = hash(rep(x), h)
convert(::Type{BigInt}, x::ZZp) = BigInt(rep(x))
Base.BigInt(x::ZZp) = rep(x)

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
    A = Array{ZZp}(undef, dims)
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
