# GF(2) and bit-packed binary polynomials --------------------------------------

export GF2, GF2X, bits, weight, NumBits, NumBytes
export random_GF2X, GF2XFactorization, BuildIrredGF2X
export to_modpoly, from_modpoly

"""The two-element field. Addition is XOR and multiplication is AND."""
struct GF2 <: Number
    bit::Bool
end
GF2(a::Integer) = GF2(isodd(mod(a, 2)))
GF2(a::GF2) = a

zero(::Type{GF2}) = GF2(false)
one(::Type{GF2}) = GF2(true)
oneunit(::Type{GF2}) = GF2(true)
zero(a::GF2) = GF2(false)
one(a::GF2) = GF2(true)
oneunit(a::GF2) = GF2(true)
iszero(a::GF2) = !a.bit
isone(a::GF2) = a.bit
rep(a::GF2) = Int(a.bit)
parent(::GF2) = GF2
characteristic(::Type{GF2}) = 2
characteristic(::GF2) = 2
modulus(::Type{GF2}) = 2
modulus(::GF2) = 2
_field_order(::Type{GF2}) = BigInt(2)

+(a::GF2, b::GF2) = GF2(xor(a.bit, b.bit))
-(a::GF2, b::GF2) = a + b
-(a::GF2) = a
*(a::GF2, b::GF2) = GF2(a.bit & b.bit)
inv(a::GF2) = iszero(a) ? throw(DomainError(a, "zero is not invertible")) : a
/(a::GF2, b::GF2) = a * inv(b)
^(a::GF2, e::Integer) = e < 0 ? inv(a) : (iszero(e) ? one(GF2) : a)
power(a::GF2, e::Integer) = a^e
sqr(a::GF2) = a

for op in (:+, :-, :*)
    @eval begin
        ($op)(a::GF2, b::Integer) = ($op)(a, GF2(b))
        ($op)(a::Integer, b::GF2) = ($op)(GF2(a), b)
    end
end
/(a::GF2, b::Integer) = a / GF2(b)
/(a::Integer, b::GF2) = GF2(a) / b
==(a::GF2, b::Integer) = rep(a) == mod(b, 2)
==(a::Integer, b::GF2) = b == a
isequal(a::GF2, b::Integer) = a == b
isequal(a::Integer, b::GF2) = b == a
hash(a::GF2, h::UInt) = hash(rep(a), h)
convert(::Type{Int}, a::GF2) = rep(a)
Base.Int(a::GF2) = rep(a)
show(io::IO, a::GF2) = print(io, "GF2(", rep(a), ")")
Random.rand(rng::AbstractRNG, ::Random.SamplerType{GF2}) = GF2(rand(rng, Bool))

"""A binary polynomial stored as nonnegative bits; bit `i` is coefficient `x^i`."""
struct GF2X
    _bits::BigInt
    function GF2X(b::BigInt)
        b >= 0 || throw(DomainError(b, "the bit representation must be nonnegative"))
        new(b)
    end
end
GF2X() = GF2X(BigInt(0))
GF2X(a::Integer) = GF2X(BigInt(a))
GF2X(a::GF2) = GF2X(BigInt(rep(a)))
function GF2X(c::AbstractVector)
    b = BigInt(0)
    for (i, a) in enumerate(c)
        !iszero(GF2(a)) && (b |= BigInt(1) << (i - 1))
    end
    GF2X(b)
end

bits(f::GF2X) = f._bits
field(::GF2X) = GF2
parent(::GF2X) = GF2
degree(f::GF2X) = iszero(f._bits) ? -1 : ndigits(f._bits, base = 2) - 1
deg(f::GF2X) = degree(f)
coeff(f::GF2X, i::Integer) = i < 0 ? throw(DomainError(i, "negative coefficient index")) : GF2(Int((f._bits >> i) & 1))
LeadCoeff(f::GF2X) = iszero(f) ? GF2(0) : GF2(1)
ConstTerm(f::GF2X) = coeff(f, 0)
coefficients(f::GF2X) = [coeff(f, i) for i in 0:degree(f)]
weight(f::GF2X) = sum(!iszero(coeff(f, i)) for i in 0:degree(f))
NumBits(f::GF2X) = degree(f) + 1
NumBytes(f::GF2X) = cld(NumBits(f), 8)
zero(::GF2X) = GF2X()
one(::GF2X) = GF2X(1)
zero(::Type{GF2X}) = GF2X()
one(::Type{GF2X}) = GF2X(1)
iszero(f::GF2X) = iszero(f._bits)
isone(f::GF2X) = f._bits == 1
==(f::GF2X, g::GF2X) = f._bits == g._bits
isequal(f::GF2X, g::GF2X) = f == g
hash(f::GF2X, h::UInt) = hash(f._bits, h)
monic(f::GF2X) = f
ismonic(f::GF2X) = !iszero(f)
show(io::IO, f::GF2X) = print(io, "GF2X(0x", string(f._bits, base = 16), ")")

+(f::GF2X, g::GF2X) = GF2X(xor(f._bits, g._bits))
-(f::GF2X, g::GF2X) = f + g
-(f::GF2X) = f

function *(f::GF2X, g::GF2X)
    a = f._bits
    b = g._bits
    result = BigInt(0)
    while b > 0
        isodd(b) && (result = xor(result, a))
        a <<= 1
        b >>= 1
    end
    GF2X(result)
end

+(f::GF2X, c::Number) = f + GF2X(GF2(c))
+(c::Number, f::GF2X) = f + c
-(f::GF2X, c::Number) = f + c
-(c::Number, f::GF2X) = f + c
*(f::GF2X, c::Number) = iszero(GF2(c)) ? GF2X() : f
*(c::Number, f::GF2X) = f * c

LeftShift(f::GF2X, n::Integer) = n >= 0 ? GF2X(f._bits << n) : throw(DomainError(n, "negative shift"))
RightShift(f::GF2X, n::Integer) = n >= 0 ? GF2X(f._bits >> n) : throw(DomainError(n, "negative shift"))
<<(f::GF2X, n::Integer) = LeftShift(f, n)
>>(f::GF2X, n::Integer) = RightShift(f, n)

function ^(f::GF2X, e::Integer)
    e >= 0 || throw(DomainError(e, "ordinary polynomial powers must be nonnegative"))
    result = one(f)
    base = f
    k = BigInt(e)
    while k > 0
        isodd(k) && (result *= base)
        k >>= 1
        k > 0 && (base *= base)
    end
    result
end

function DivRem(f::GF2X, g::GF2X)
    iszero(g) && throw(DivideError())
    r = f._bits
    q = BigInt(0)
    dg = degree(g)
    while r > 0
        dr = iszero(r) ? -1 : ndigits(r, base = 2) - 1
        dr >= dg || break
        shift = dr - dg
        q |= BigInt(1) << shift
        r = xor(r, g._bits << shift)
    end
    GF2X(q), GF2X(r)
end
div(f::GF2X, g::GF2X) = DivRem(f, g)[1]
rem(f::GF2X, g::GF2X) = DivRem(f, g)[2]
/(f::GF2X, g::GF2X) = div(f, g)
/(f::GF2X, c::Number) = iszero(GF2(c)) ? throw(DivideError()) : f

function gcd(f::GF2X, g::GF2X)
    a, b = f, g
    while !iszero(b)
        a, b = b, rem(a, b)
    end
    a
end
GCD(f::GF2X, g::GF2X) = gcd(f, g)

function XGCD(f::GF2X, g::GF2X)
    old_r, r = f, g
    old_s, s = one(f), zero(f)
    old_t, t = zero(f), one(f)
    while !iszero(r)
        q = div(old_r, r)
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    end
    old_r, old_s, old_t
end

function derivative(f::GF2X)
    result = BigInt(0)
    for i in 1:2:degree(f)
        !iszero(coeff(f, i)) && (result |= BigInt(1) << (i - 1))
    end
    GF2X(result)
end
diff(f::GF2X) = derivative(f)
evaluate(f::GF2X, x) = iszero(GF2(x)) ? ConstTerm(f) : GF2(weight(f) % 2)
compose(f::GF2X, g::GF2X) = begin
    result = zero(f)
    for i in degree(f):-1:0
        result = result * g + coeff(f, i)
    end
    result
end

function reverse(f::GF2X)
    iszero(f) && return f
    result = BigInt(0)
    for i in 0:degree(f)
        !iszero(coeff(f, i)) && (result |= BigInt(1) << (degree(f) - i))
    end
    GF2X(result)
end
truncate(f::GF2X, n::Integer) = begin
    n >= 0 || throw(DomainError(n, "negative truncation"))
    n == 0 ? GF2X() : GF2X(f._bits & ((BigInt(1) << n) - 1))
end
MulTrunc(f::GF2X, g::GF2X, n::Integer) = truncate(f * g, n)
MulMod(f::GF2X, g::GF2X, m::GF2X) = rem(f * g, m)

function InvMod(f::GF2X, m::GF2X)
    d, s, _ = XGCD(f, m)
    isone(d) || throw(DomainError(f, "polynomial inverse undefined"))
    rem(s, m)
end
function InvModStatus(f::GF2X, m::GF2X)
    d, s, _ = XGCD(f, m)
    isone(d) ? (true, rem(s, m)) : (false, zero(f))
end

function PowerMod(f::GF2X, e::Integer, m::GF2X)
    base = rem(f, m)
    k = BigInt(e)
    if k < 0
        base = InvMod(base, m)
        k = -k
    end
    result = rem(one(f), m)
    while k > 0
        isodd(k) && (result = MulMod(result, base, m))
        k >>= 1
        k > 0 && (base = MulMod(base, base, m))
    end
    result
end

function InvTrunc(f::GF2X, n::Integer)
    n >= 0 || throw(DomainError(n, "negative truncation"))
    n == 0 && return GF2X()
    iszero(ConstTerm(f)) && throw(DomainError(f, "constant term must be one"))
    g = GF2X(1)
    m = 1
    while m < n
        m2 = min(2m, Int(n))
        # In characteristic two, 2-f*g equals f*g, so Newton doubling
        # becomes g' = f*g^2 rather than the char-0-looking 1+f*g form.
        g = truncate(truncate(f, m2) * g * g, m2)
        m = m2
    end
    g
end

random_GF2X(n::Integer; rng::AbstractRNG = Random.default_rng()) = begin
    n >= 0 || throw(DomainError(n, "negative coefficient count"))
    GF2X(RandomBits(n; rng = rng))
end

to_modpoly(f::GF2X) = poly(GF2, coefficients(f))
from_modpoly(f::ModPoly) = field(f) == GF2 ? GF2X(coefficients(f)) : throw(ArgumentError("expected coefficients over GF2"))

struct GF2XFactorization
    unit::GF2
    factors::Vector{Tuple{GF2X,Int}}
end

function factor(f::GF2X; rng::AbstractRNG = Random.default_rng())
    iszero(f) && throw(DomainError(f, "zero has no factorization"))
    fac = factor(to_modpoly(f); rng = rng)
    GF2XFactorization(fac.unit, [(from_modpoly(g), m) for (g, m) in fac.factors])
end
function reconstruct(fac::GF2XFactorization)
    result = GF2X(fac.unit)
    for (g, m) in fac.factors
        result *= g^m
    end
    result
end
is_irreducible(f::GF2X) = is_irreducible(to_modpoly(f))
BuildIrredGF2X(n::Integer; rng::AbstractRNG = Random.default_rng()) = from_modpoly(BuildIrred(GF2, n; rng = rng))
FindRoots(f::GF2X; rng::AbstractRNG = Random.default_rng()) = FindRoots(to_modpoly(f); rng = rng)
FindRoot(f::GF2X; rng::AbstractRNG = Random.default_rng()) = FindRoot(to_modpoly(f); rng = rng)
