module NTL

using Random

import Base: +, -, *, /, ^, inv, zero, one, oneunit, iszero, isone
import Base: ==, isequal, hash, convert, promote_rule, show, parent, diff, truncate, muladd

export PrimeField, zzp, zzp_Min, zzp_Max
export rep, modulus, characteristic, info
export AddMod, SubMod, NegateMod, MulMod, SqrMod
export InvMod, InvModStatus, PowerMod
export add, sub, negate, mul, sqr, divide, power, muladd
export random_zzp, VectorRandom, VectorConv, VectorMulMod, VectorMulMod!, swap

# Internal marker used to construct an already reduced representative.
struct _UnsafeRep end
const _unsafe_rep = _UnsafeRep()

"""
    PrimeField{p} <: Number

An element of the prime field `F_p`. The modulus is a Julia type parameter,
so values with different moduli cannot be mixed accidentally. Elements are
immutable and always store a representative in `0:p-1`.

Use `F = PrimeField(p)` to validate `p` and obtain the field element type, then use
`F(a)` to construct an element. Direct `PrimeField{p}(a)` construction is also
supported and performs the same primality check. Iterating `PrimeField{p}` yields
all field elements in canonical order, from `0` through `p-1`.
"""
struct PrimeField{P} <: Number
    _rep::Int
    function PrimeField{P}(r::Int, ::_UnsafeRep) where {P}
        new{P}(r)
    end
end

# Legacy NTL spelling; new code should use PrimeField.
const zzp = PrimeField

@inline _raw(::Val{P}, r::Int) where {P} = PrimeField{P}(r, _unsafe_rep)

const _MULMOD_DIRECT_MAX = isqrt(typemax(Int))

# Use a machine-word product whenever it cannot overflow. Wider UInt128 is used
# only for larger moduli, where the product is exact as well.
@inline function _mulmod_raw(a::Int, b::Int, n::Int)
    n <= _MULMOD_DIRECT_MAX ? (a * b) % n : Int(mod(UInt128(a) * UInt128(b), UInt128(n)))
end

@inline _muladd_raw(a::Int, b::Int, c::Int, n::Int) = _addmod_raw(_mulmod_raw(a, b, n), c, n)

# Overflow-free addition. The branch guarantees that a+b < n when it is used.
@inline function _addmod_raw(a::Int, b::Int, n::Int)
    t = n - b
    a >= t ? a - t : a + b
end

# Overflow-free subtraction. In the second branch, a+(n-b) < n.
@inline function _submod_raw(a::Int, b::Int, n::Int)
    a >= b ? a - b : a + (n - b)
end

@inline _negatemod_raw(a::Int, n::Int) = iszero(a) ? 0 : n - a

function _powermod_uint(a::Int, e::UInt, n::Int)
    x = mod(1, n)
    y = a
    k = e
    while !iszero(k)
        if isodd(k)
            x = _mulmod_raw(x, y, n)
        end
        k >>= 1
        if !iszero(k)
            y = _mulmod_raw(y, y, n)
        end
    end
    x
end

function _powermod_big(a::Int, e::BigInt, n::Int)
    x = mod(1, n)
    y = a
    k = e
    while !iszero(k)
        if isodd(k)
            x = _mulmod_raw(x, y, n)
        end
        k >>= 1
        if !iszero(k)
            y = _mulmod_raw(y, y, n)
        end
    end
    x
end

# Deterministic Miller-Rabin for positive Int values. These witnesses are
# deterministic below 2^64; products are evaluated in Int128.
function _is_prime(n::Int)
    n < 2 && return false
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
        n == p && return true
        n % p == 0 && return false
    end

    d = n - 1
    s = 0
    while iseven(d)
        d >>= 1
        s += 1
    end

    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022)
        a % n == 0 && continue
        x = _powermod_uint(a % n, UInt(d), n)
        (x == 1 || x == n - 1) && continue
        passed = false
        for _ in 1:(s - 1)
            x = _mulmod_raw(x, x, n)
            if x == n - 1
                passed = true
                break
            end
        end
        if x != n - 1 && !passed
            return false
        end
    end
    true
end

function _previous_prime(n::Int)
    n >= 2 || throw(DomainError(n, "no prime is below this bound"))
    k = n
    iseven(k) && (k -= 1)
    while k >= 2
        _is_prime(k) && return k
        k -= 2
    end
    2
end

"""Minimum supported prime modulus."""
const zzp_Min = 2

"""Largest prime representable by the platform's signed `Int` type."""
const zzp_Max = _previous_prime(typemax(Int))

function _require_prime(p::Int)
    p >= zzp_Min || throw(DomainError(p, "the modulus must be at least 2"))
    _is_prime(p) || throw(ArgumentError("the modulus $p is not prime"))
    p
end

@generated function _checked_modulus(::Val{P}) where {P}
    P isa Int || return :(throw(ArgumentError("the modulus type parameter must be an Int")))
    P < zzp_Min && return :(throw(DomainError($P, "the modulus must be at least 2")))
    _is_prime(P) || return :(throw(ArgumentError("the modulus is not prime")))
    P
end

@inline _from_integer(::Val{P}, a::Integer) where {P} = _raw(Val{P}(), Int(mod(a, P)))

"""
    PrimeField(p::Integer)

Validate the prime `p` and return the element type `PrimeField{p}`.
"""
function PrimeField(p::Integer)
    2 <= p <= typemax(Int) || throw(DomainError(p, "the modulus must fit in Int and be at least 2"))
    P = Int(p)
    _require_prime(P)
    PrimeField{P}
end

PrimeField{P}() where {P} = (_checked_modulus(Val{P}()); _raw(Val{P}(), 0))
PrimeField{P}(a::Integer) where {P} = (_checked_modulus(Val{P}()); _from_integer(Val{P}(), a))
PrimeField{P}(a::PrimeField{P}) where {P} = a
function PrimeField{P}(a::PrimeField) where {P}
    throw(ArgumentError("cannot convert an element from modulus $(modulus(a)) to modulus $P"))
end

"""Return the canonical integer representative in `0:p-1`."""
rep(a::PrimeField) = a._rep

"""Return the prime modulus of a field type or element."""
modulus(::Type{PrimeField{P}}) where {P} = P
modulus(::PrimeField{P}) where {P} = P

"""Return the field characteristic; equivalent to `modulus`."""
characteristic(F::Type{PrimeField{P}}) where {P} = modulus(F)
characteristic(a::PrimeField) = modulus(a)

"""
    for a in PrimeField{p} ... end

Iterate over all `p` field elements in canonical order. `collect(PrimeField{p})`
returns the same sequence as a vector; avoid collecting very large fields.
"""
function Base.iterate(F::Type{PrimeField{P}}, state::Int = 0) where {P}
    state == 0 && _checked_modulus(Val{P}())
    state < P || return nothing
    (_raw(Val{P}(), state), state + 1)
end
Base.length(::Type{PrimeField{P}}) where {P} = (_checked_modulus(Val{P}()); P)
Base.eltype(::Type{PrimeField{P}}) where {P} = PrimeField{P}
Base.IteratorSize(::Type{PrimeField{P}}) where {P} = Base.HasLength()
Base.IteratorEltype(::Type{PrimeField{P}}) where {P} = Base.HasEltype()
Base.first(F::Type{PrimeField{P}}) where {P} = zero(F)
Base.last(F::Type{PrimeField{P}}) where {P} = (_checked_modulus(Val{P}()); _raw(Val{P}(), P - 1))

"""Return summary information about the implementation or a specific field."""
info(::Type{PrimeField}) = (
    description = "Single-word prime-field arithmetic inspired by NTL lzz_p/zz_p",
    minimum_modulus = zzp_Min,
    maximum_modulus = zzp_Max,
    representation = "canonical signed-Int residue",
)
function info(F::Type{PrimeField{P}}) where {P}
    (modulus = P, zero = zero(F), one = one(F), storage = Int)
end
info(a::PrimeField{P}) where {P} = (representative = rep(a), modulus = P, field = PrimeField{P})

zero(F::Type{PrimeField{P}}) where {P} = F(0)
zero(a::PrimeField) = zero(typeof(a))
one(F::Type{PrimeField{P}}) where {P} = F(1)
one(a::PrimeField) = one(typeof(a))
oneunit(F::Type{PrimeField{P}}) where {P} = one(F)
oneunit(a::PrimeField) = one(a)

iszero(a::PrimeField) = iszero(rep(a))
isone(a::PrimeField) = isone(rep(a))

convert(F::Type{PrimeField{P}}, a::Integer) where {P} = F(a)
convert(::Type{PrimeField{P}}, a::PrimeField{P}) where {P} = a
function convert(::Type{PrimeField{P}}, a::PrimeField) where {P}
    throw(ArgumentError("cannot convert between different prime fields"))
end
convert(::Type{T}, a::PrimeField) where {T<:Integer} = convert(T, rep(a))
Base.Int(a::PrimeField) = rep(a)
Base.BigInt(a::PrimeField) = BigInt(rep(a))
promote_rule(::Type{PrimeField{P}}, ::Type{T}) where {P,T<:Integer} = PrimeField{P}

show(io::IO, a::PrimeField{P}) where {P} = print(io, "PrimeField{", P, "}(", rep(a), ")")
show(io::IO, ::Type{PrimeField{P}}) where {P} = print(io, "PrimeField{", P, "}")

# Low-level modular arithmetic -------------------------------------------------

function _check_positive_modulus(n::Int)
    n > 0 || throw(DomainError(n, "the modulus must be positive"))
    n
end

function _check_rep(a::Int, n::Int)
    0 <= a < n || throw(DomainError(a, "expected a representative in 0:n-1"))
    a
end

"""Return `(a+b) mod n` for reduced `Int` representatives `a` and `b`."""
function AddMod(a::Int, b::Int, n::Int)
    _check_positive_modulus(n)
    _check_rep(a, n); _check_rep(b, n)
    _addmod_raw(a, b, n)
end

"""Return `(a-b) mod n` for reduced `Int` representatives `a` and `b`."""
function SubMod(a::Int, b::Int, n::Int)
    _check_positive_modulus(n)
    _check_rep(a, n); _check_rep(b, n)
    _submod_raw(a, b, n)
end

"""Return `(-a) mod n` for a reduced `Int` representative `a`."""
function NegateMod(a::Int, n::Int)
    _check_positive_modulus(n)
    _check_rep(a, n)
    _negatemod_raw(a, n)
end

"""Return `(a*b) mod n` using an intermediate `Int128`."""
function MulMod(a::Int, b::Int, n::Int)
    _check_positive_modulus(n)
    _check_rep(a, n); _check_rep(b, n)
    _mulmod_raw(a, b, n)
end

"""Return `(a^2) mod n` using an intermediate `Int128`."""
SqrMod(a::Int, n::Int) = MulMod(a, a, n)

"""
    InvModStatus(a::Int, n::Int) -> (Bool, Int)

Attempt to compute an inverse modulo `n`. Returns `(true, inverse)` when
`gcd(a,n)==1`, and `(false, 0)` otherwise.
"""
function InvModStatus(a::Int, n::Int)
    n > 1 || throw(DomainError(n, "an inverse requires modulus greater than 1"))
    r = mod(a, n)
    g, u, _ = gcdx(r, n)
    g == 1 ? (true, mod(u, n)) : (false, 0)
end

@inline function _invmod_reduced(a::Int, n::Int)
    g, u, _ = gcdx(a, n)
    g == 1 || throw(DomainError(a, "inverse undefined: gcd(a,n) is not 1"))
    u < 0 ? u + n : u
end

"""Return the multiplicative inverse of `a` modulo `n`, or throw `DomainError`."""
function InvMod(a::Int, n::Int)
    n > 1 || throw(DomainError(n, "an inverse requires modulus greater than 1"))
    r = mod(a, n)
    _invmod_reduced(r, n)
end

"""Return `(a^e) mod n` by binary square-and-multiply; negative exponents invert `a`."""
function PowerMod(a::Int, e::Int, n::Int)
    _check_positive_modulus(n)
    base = mod(a, n)
    if e < 0
        base = InvMod(base, n)
        k = UInt(-(e + 1)) + UInt(1)
        return _powermod_uint(base, k, n)
    end
    _powermod_uint(base, UInt(e), n)
end

function PowerMod(a::Int, e::Integer, n::Int)
    _check_positive_modulus(n)
    base = mod(a, n)
    k = BigInt(e)
    if k < 0
        base = InvMod(base, n)
        k = -k
    end
    _powermod_big(base, k, n)
end

# Field arithmetic -------------------------------------------------------------

+(x::PrimeField{P}, y::PrimeField{P}) where {P} = _raw(Val{P}(), _addmod_raw(rep(x), rep(y), P))
+(x::PrimeField, y::PrimeField) = throw(ArgumentError("cannot add elements with different moduli"))
+(x::PrimeField{P}, y::Integer) where {P} = x + _from_integer(Val{P}(), y)
+(x::Integer, y::PrimeField{P}) where {P} = _from_integer(Val{P}(), x) + y

-(x::PrimeField{P}, y::PrimeField{P}) where {P} = _raw(Val{P}(), _submod_raw(rep(x), rep(y), P))
-(x::PrimeField, y::PrimeField) = throw(ArgumentError("cannot subtract elements with different moduli"))
-(x::PrimeField{P}, y::Integer) where {P} = x - _from_integer(Val{P}(), y)
-(x::Integer, y::PrimeField{P}) where {P} = _from_integer(Val{P}(), x) - y
-(x::PrimeField{P}) where {P} = _raw(Val{P}(), _negatemod_raw(rep(x), P))

*(x::PrimeField{P}, y::PrimeField{P}) where {P} = _raw(Val{P}(), _mulmod_raw(rep(x), rep(y), P))
*(x::PrimeField, y::PrimeField) = throw(ArgumentError("cannot multiply elements with different moduli"))
*(x::PrimeField{P}, y::Integer) where {P} = x * _from_integer(Val{P}(), y)
*(x::Integer, y::PrimeField{P}) where {P} = _from_integer(Val{P}(), x) * y

@inline muladd(x::PrimeField{P}, y::PrimeField{P}, z::PrimeField{P}) where {P} = _raw(Val{P}(), _muladd_raw(rep(x), rep(y), rep(z), P))

inv(x::PrimeField{P}) where {P} = _raw(Val{P}(), _invmod_reduced(rep(x), P))
/(x::PrimeField{P}, y::PrimeField{P}) where {P} = x * inv(y)
/(x::PrimeField, y::PrimeField) = throw(ArgumentError("cannot divide elements with different moduli"))
/(x::PrimeField{P}, y::Integer) where {P} = x / _from_integer(Val{P}(), y)
/(x::Integer, y::PrimeField{P}) where {P} = _from_integer(Val{P}(), x) / y

function power(x::PrimeField{P}, e::Integer) where {P}
    e >= 0 && return Base.power_by_squaring(x, e)
    e == typemin(Int) && return _raw(Val{P}(), PowerMod(rep(x), e, P))
    Base.power_by_squaring(inv(x), -e)
end
^(x::PrimeField, e::Integer) = power(x, e)

==(x::PrimeField{P}, y::PrimeField{P}) where {P} = rep(x) == rep(y)
==(::PrimeField, ::PrimeField) = false
==(x::PrimeField{P}, y::Integer) where {P} = rep(x) == mod(y, P)
==(x::Integer, y::PrimeField) = y == x
isequal(x::PrimeField{P}, y::PrimeField{P}) where {P} = rep(x) == rep(y)
isequal(::PrimeField, ::PrimeField) = false
isequal(x::PrimeField, y::Integer) = x == y
isequal(x::Integer, y::PrimeField) = y == x
hash(x::PrimeField, h::UInt) = hash(rep(x), h)

# Procedural forms mirror the algorithms while returning fresh immutable values.
add(x, y) = x + y
sub(x, y) = x - y
negate(x) = -x
mul(x, y) = x * y
sqr(x::PrimeField{P}) where {P} = _raw(Val{P}(), _mulmod_raw(rep(x), rep(x), P))
divide(x, y) = x / y

"""Return a uniform random field element; equivalent to `rand(F)`."""
random_zzp(F::Type{PrimeField{P}}; rng::AbstractRNG = Random.default_rng()) where {P} = rand(rng, F)

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{PrimeField{P}}) where {P}
    _checked_modulus(Val{P}())
    _raw(Val{P}(), rand(rng, 0:(P - 1)))
end

"""Return a vector of `k` random elements of `F`."""
function VectorRandom(F::Type{PrimeField{P}}, k::Integer; rng::AbstractRNG = Random.default_rng()) where {P}
    k >= 0 || throw(DomainError(k, "the vector length must be nonnegative"))
    rand(rng, F, Int(k))
end

"""Convert an integer vector to canonical elements of `F`."""
VectorConv(F::Type{PrimeField{P}}, a::AbstractVector{<:Integer}) where {P} = [F(x) for x in a]

"""Return `a*b mod n` componentwise for reduced integer representatives."""
function VectorMulMod(a::AbstractVector{Int}, b::Int, n::Int)
    _check_positive_modulus(n)
    _check_rep(b, n)
    out = Vector{Int}(undef, length(a))
    for i in eachindex(a)
        _check_rep(a[i], n)
        out[i] = _mulmod_raw(a[i], b, n)
    end
    out
end

"""Write `a*b mod n` componentwise into `x`; aliasing `x === a` is safe."""
function VectorMulMod!(x::AbstractVector{Int}, a::AbstractVector{Int}, b::Int, n::Int)
    length(x) == length(a) || throw(DimensionMismatch("input and output lengths differ"))
    _check_positive_modulus(n)
    _check_rep(b, n)
    for i in eachindex(x, a)
        _check_rep(a[i], n)
        x[i] = _mulmod_raw(a[i], b, n)
    end
    x
end

"""Return the pair `(y,x)`; Julia's immutable field elements are not swapped in place."""
swap(x::PrimeField{P}, y::PrimeField{P}) where {P} = (y, x)

include("zzp_big.jl")
include("numbertheory.jl")
include("polynomial.jl")
include("binary.jl")
include("extension.jl")
include("zech.jl")
include("matrix.jl")
include("lattice.jl")
include("lll_fast.jl")
include("transform.jl")
include("fft.jl")
include("polynomial_mul.jl")

end # module NTL
