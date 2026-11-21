# Arbitrary-precision number theory -------------------------------------------

export ProbPrime, MillerWitness, NextPrime, PrevPrime, GenPrime
export RandomBits, RandomLen, RandomBnd, FactorInt, EulerPhi
export CRT, SqrRootMod, DiscreteLog, PrimitiveRoot

function _strong_mr_pass(n::BigInt, a::BigInt, d::BigInt, s::Int)
    a = mod(a, n)
    iszero(a) && return true
    x = powermod(a, d, n)
    (x == 1 || x == n - 1) && return true
    for _ in 1:(s - 1)
        x = mod(x * x, n)
        x == n - 1 && return true
    end
    false
end

"""Miller–Rabin probable-prime test; deterministic in the signed-Int range."""
function ProbPrime(n::Integer; rounds::Integer = 16)
    rounds > 0 || throw(DomainError(rounds, "rounds must be positive"))
    n < 2 && return false
    if n <= typemax(Int)
        return _is_prime(Int(n))
    end
    N = BigInt(n)
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small
        N == p && return true
        N % p == 0 && return false
    end
    d = N - 1
    s = 0
    while iseven(d)
        d >>= 1
        s += 1
    end
    for i in 1:Int(rounds)
        a = i <= length(small) ? BigInt(small[i]) : BigInt(i + 35)
        _strong_mr_pass(N, a, d, s) || return false
    end
    true
end

"""Return `true` when `w` is a Miller witness for the compositeness of odd `n`."""
function MillerWitness(n::Integer, w::Integer)
    n > 2 && isodd(n) || throw(DomainError(n, "n must be an odd integer greater than 2"))
    N = BigInt(n)
    d = N - 1
    s = 0
    while iseven(d)
        d >>= 1
        s += 1
    end
    !_strong_mr_pass(N, BigInt(w), d, s)
end

"""Return the least probable prime greater than or equal to `m`."""
function NextPrime(m::Integer; rounds::Integer = 16)
    m <= 2 && return BigInt(2)
    n = BigInt(m)
    iseven(n) && (n += 1)
    while !ProbPrime(n; rounds = rounds)
        n += 2
    end
    n
end

"""Return the greatest probable prime less than or equal to `m`."""
function PrevPrime(m::Integer; rounds::Integer = 16)
    m >= 2 || throw(DomainError(m, "no prime is below this bound"))
    m < 3 && return BigInt(2)
    n = BigInt(m)
    iseven(n) && (n -= 1)
    while n >= 3 && !ProbPrime(n; rounds = rounds)
        n -= 2
    end
    n >= 3 ? n : BigInt(2)
end

"""Return a uniformly generated nonnegative BigInt below `n`."""
RandomBnd(n::Integer; rng::AbstractRNG = Random.default_rng()) = _random_big_below(rng, BigInt(n))

"""Return a random nonnegative integer with at most `bits` random bits."""
function RandomBits(bits::Integer; rng::AbstractRNG = Random.default_rng())
    bits >= 0 || throw(DomainError(bits, "the bit count must be nonnegative"))
    bits == 0 && return BigInt(0)
    bound = BigInt(1) << Int(bits)
    _random_big_below(rng, bound)
end

"""Return a random integer with exactly `bits` bits and its top bit set."""
function RandomLen(bits::Integer; rng::AbstractRNG = Random.default_rng())
    bits > 0 || throw(DomainError(bits, "the bit count must be positive"))
    bits == 1 && return BigInt(1)
    (BigInt(1) << Int(bits - 1)) | RandomBits(bits - 1; rng = rng)
end

"""Generate a probable prime with exactly `bits` bits."""
function GenPrime(bits::Integer; rounds::Integer = 16, rng::AbstractRNG = Random.default_rng())
    bits >= 2 || throw(DomainError(bits, "a prime needs at least two bits"))
    while true
        n = RandomLen(bits; rng = rng) | 1
        ProbPrime(n; rounds = rounds) && return n
    end
end

# Brent's batch-gcd variant of Pollard rho.
function _pollard_rho(n::BigInt)
    n % 2 == 0 && return BigInt(2)
    c = BigInt(1)
    while true
        y = BigInt(2)
        m = BigInt(128)
        g = BigInt(1)
        r = BigInt(1)
        q = BigInt(1)
        x = BigInt(0)
        ys = BigInt(0)
        while g == 1
            x = y
            for _ in 1:r
                y = mod(y * y + c, n)
            end
            k = BigInt(0)
            while k < r && g == 1
                ys = y
                for _ in 1:min(m, r - k)
                    y = mod(y * y + c, n)
                    q = mod(q * abs(x - y), n)
                end
                g = gcd(q, n)
                k += m
            end
            r <<= 1
        end
        if g == n
            while true
                ys = mod(ys * ys + c, n)
                g = gcd(abs(x - ys), n)
                g > 1 && break
            end
        end
        g != n && return g
        c += 1
    end
end

function _factor_recursive(n::BigInt, out::Vector{BigInt})
    n == 1 && return
    if ProbPrime(n)
        push!(out, n)
        return
    end
    d = _pollard_rho(n)
    _factor_recursive(d, out)
    _factor_recursive(div(n, d), out)
end

"""
    FactorInt(n) -> Vector{Tuple{BigInt,Int}}

Factor a nonzero integer into `(prime, exponent)` pairs. A negative input starts
with `(-1,1)`. The implementation combines Miller–Rabin and Pollard–Brent rho.
"""
function FactorInt(n::Integer)
    iszero(n) && throw(DomainError(n, "zero has no prime factorization"))
    N = BigInt(n)
    negative = N < 0
    N = abs(N)
    pieces = BigInt[]
    N > 1 && _factor_recursive(N, pieces)
    sort!(pieces)
    out = Tuple{BigInt,Int}[]
    negative && push!(out, (BigInt(-1), 1))
    for p in pieces
        if !isempty(out) && out[end][1] == p
            out[end] = (p, out[end][2] + 1)
        else
            push!(out, (p, 1))
        end
    end
    out
end

"""Compute Euler's totient from `FactorInt`."""
function EulerPhi(n::Integer)
    n > 0 || throw(DomainError(n, "EulerPhi expects a positive integer"))
    n == 1 && return BigInt(1)
    result = BigInt(n)
    for (p, _) in FactorInt(n)
        result = div(result, p) * (p - 1)
    end
    result
end

"""Generalized Chinese remaindering for pairwise coprime positive moduli."""
function CRT(residues::AbstractVector{<:Integer}, moduli::AbstractVector{<:Integer})
    length(residues) == length(moduli) || throw(DimensionMismatch("residue and modulus lengths differ"))
    isempty(moduli) && throw(ArgumentError("CRT needs at least one congruence"))
    x = BigInt(0)
    M = BigInt(1)
    for (r, m) in zip(residues, moduli)
        m > 0 || throw(DomainError(m, "CRT moduli must be positive"))
        gcd(M, m) == 1 || throw(ArgumentError("CRT moduli must be pairwise coprime"))
        step = mod((BigInt(r) - x) * invmod(M, m), m)
        x = mod(x + M * step, M * m)
        M *= m
    end
    x, M
end

"""Tonelli–Shanks square root modulo a prime; return `nothing` for a nonresidue."""
function SqrRootMod(a::Integer, p::Integer)
    P = BigInt(p)
    ProbPrime(P) || throw(ArgumentError("the modulus must be prime"))
    A = mod(BigInt(a), P)
    P == 2 && return A
    iszero(A) && return BigInt(0)
    powermod(A, div(P - 1, 2), P) == 1 || return nothing
    P % 4 == 3 && return powermod(A, div(P + 1, 4), P)

    q = P - 1
    s = 0
    while iseven(q)
        q >>= 1
        s += 1
    end
    z = BigInt(2)
    while powermod(z, div(P - 1, 2), P) != P - 1
        z += 1
    end
    c = powermod(z, q, P)
    x = powermod(A, div(q + 1, 2), P)
    t = powermod(A, q, P)
    m = s
    while t != 1
        i = 1
        t2 = mod(t * t, P)
        while i < m && t2 != 1
            t2 = mod(t2 * t2, P)
            i += 1
        end
        i < m || error("Tonelli–Shanks failed to find a square-root order")
        b = powermod(c, BigInt(1) << (m - i - 1), P)
        b2 = mod(b * b, P)
        x = mod(x * b, P)
        t = mod(t * b2, P)
        c = b2
        m = i
    end
    x
end

"""Baby-step/giant-step discrete logarithm `x` with `a^x = b mod p`."""
function DiscreteLog(a::Integer, b::Integer, p::Integer)
    P = BigInt(p)
    ProbPrime(P) || throw(ArgumentError("the modulus must be prime"))
    A = mod(BigInt(a), P)
    B = mod(BigInt(b), P)
    gcd(A, P) == 1 || throw(DomainError(a, "the base must be a unit"))
    P == 2 && return B == 1 ? BigInt(0) : nothing
    mbig = isqrt(P - 1) + 1
    mbig <= typemax(Int) || throw(ArgumentError("the required table is too large"))
    m = Int(mbig)
    table = Dict{BigInt,Int}()
    e = BigInt(1)
    for j in 0:(m - 1)
        get!(table, e, j)
        e = mod(e * A, P)
    end
    factor = invmod(powermod(A, m, P), P)
    gamma = B
    for i in 0:m
        if haskey(table, gamma)
            x = BigInt(i) * m + table[gamma]
            x < P - 1 && powermod(A, x, P) == B && return x
        end
        gamma = mod(gamma * factor, P)
    end
    nothing
end

"""Find a generator of the multiplicative group modulo an odd or even prime."""
function PrimitiveRoot(p::Integer)
    P = BigInt(p)
    ProbPrime(P) || throw(ArgumentError("the modulus must be prime"))
    P == 2 && return BigInt(1)
    factors = FactorInt(P - 1)
    g = BigInt(2)
    while g < P
        if all(powermod(g, div(P - 1, q), P) != 1 for (q, _) in factors)
            return g
        end
        g += 1
    end
    error("no primitive root found")
end
