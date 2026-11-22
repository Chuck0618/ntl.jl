# Arbitrary-precision modular integers ----------------------------------------

export BigModRing, ModRing, ModInt, ZZPField, ZZp, ZZ_p, ZZP, parent

"""A reusable parent object for arithmetic modulo an arbitrary positive BigInt."""
struct BigModRing
    _modulus::BigInt
    function BigModRing(n::Integer)
        n > 0 || throw(DomainError(n, "the modulus must be positive"))
        new(BigInt(n))
    end
end

# NTL-compatible spellings; new code should use BigModRing and ModInt.
const ModRing = BigModRing
const ZZPField = BigModRing
const ZZP = BigModRing

"""An arbitrary-precision modular integer. Use `R = BigModRing(n); R(a)` or `ModInt(R, a)`."""
struct ModInt <: Number
    _parent::BigModRing
    _rep::BigInt
    function ModInt(F::BigModRing, a::Integer)
        new(F, mod(BigInt(a), modulus(F)))
    end
end

const ZZp = ModInt
const ZZ_p = ModInt

(F::BigModRing)(a::Integer) = ModInt(F, a)
(F::BigModRing)(a::ModInt) = (parent(a) == F || throw(ArgumentError("modulus mismatch")); a)

modulus(F::BigModRing) = F._modulus
modulus(a::ModInt) = modulus(parent(a))
characteristic(F::BigModRing) = modulus(F)
characteristic(a::ModInt) = modulus(a)
parent(a::ModInt) = a._parent
parent(::PrimeField{P}) where {P} = PrimeField{P}
rep(a::ModInt) = a._rep

==(F::BigModRing, G::BigModRing) = modulus(F) == modulus(G)
isequal(F::BigModRing, G::BigModRing) = F == G
hash(F::BigModRing, h::UInt) = hash(modulus(F), h)

zero(F::BigModRing) = F(0)
one(F::BigModRing) = F(1)
oneunit(F::BigModRing) = F(1)
zero(a::ModInt) = zero(parent(a))
one(a::ModInt) = one(parent(a))
oneunit(a::ModInt) = one(a)
iszero(a::ModInt) = iszero(rep(a))
isone(a::ModInt) = isone(rep(a))

function show(io::IO, F::BigModRing)
    print(io, "BigModRing(", modulus(F), ")")
end
function show(io::IO, a::ModInt)
    print(io, "ModInt(", modulus(a), ")(", rep(a), ")")
end

function _check_same_parent(x::ModInt, y::ModInt)
    parent(x) == parent(y) || throw(ArgumentError("cannot mix different moduli"))
    parent(x)
end

+(x::ModInt, y::ModInt) = (F = _check_same_parent(x, y); F(rep(x) + rep(y)))
-(x::ModInt, y::ModInt) = (F = _check_same_parent(x, y); F(rep(x) - rep(y)))
-(x::ModInt) = parent(x)(-rep(x))
*(x::ModInt, y::ModInt) = (F = _check_same_parent(x, y); F(rep(x) * rep(y)))

+(x::ModInt, y::Integer) = x + parent(x)(y)
+(x::Integer, y::ModInt) = parent(y)(x) + y
-(x::ModInt, y::Integer) = x - parent(x)(y)
-(x::Integer, y::ModInt) = parent(y)(x) - y
*(x::ModInt, y::Integer) = x * parent(x)(y)
*(x::Integer, y::ModInt) = parent(y)(x) * y

function inv(x::ModInt)
    n = modulus(x)
    gcd(rep(x), n) == 1 || throw(DomainError(x, "inverse undefined for a nonunit"))
    parent(x)(invmod(rep(x), n))
end
/(x::ModInt, y::ModInt) = x * inv(y)
/(x::ModInt, y::Integer) = x / parent(x)(y)
/(x::Integer, y::ModInt) = parent(y)(x) / y
sqr(x::ModInt) = x * x

function power(x::ModInt, e::Integer)
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
^(x::ModInt, e::Integer) = power(x, e)

==(x::ModInt, y::ModInt) = parent(x) == parent(y) && rep(x) == rep(y)
isequal(x::ModInt, y::ModInt) = x == y
==(x::ModInt, y::Integer) = rep(x) == mod(BigInt(y), modulus(x))
==(x::Integer, y::ModInt) = y == x
isequal(x::ModInt, y::Integer) = x == y
isequal(x::Integer, y::ModInt) = y == x
hash(x::ModInt, h::UInt) = hash(rep(x), h)
convert(::Type{BigInt}, x::ModInt) = BigInt(rep(x))
Base.BigInt(x::ModInt) = rep(x)

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

Random.rand(F::BigModRing) = Random.rand(Random.default_rng(), F)
Random.rand(rng::AbstractRNG, F::BigModRing) = F(_random_big_below(rng, modulus(F)))
Random.rand(F::BigModRing, dims::Dims) = Random.rand(Random.default_rng(), F, dims)
Random.rand(F::BigModRing, d::Integer, dims::Integer...) = Random.rand(Random.default_rng(), F, (d, dims...))
function Random.rand(rng::AbstractRNG, F::BigModRing, dims::Dims)
    A = Array{ModInt}(undef, dims)
    for i in eachindex(A)
        A[i] = rand(rng, F)
    end
    A
end
function Random.rand(rng::AbstractRNG, F::BigModRing, d::Integer, dims::Integer...)
    rand(rng, F, (d, dims...))
end
random_zzp(F::BigModRing; rng::AbstractRNG = Random.default_rng()) = rand(rng, F)
VectorRandom(F::BigModRing, k::Integer; rng::AbstractRNG = Random.default_rng()) = rand(rng, F, Int(k))
VectorConv(F::BigModRing, a::AbstractVector{<:Integer}) = [F(x) for x in a]
