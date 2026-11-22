# Dense univariate polynomials -------------------------------------------------

import Base: div, rem, gcd, %, <<, >>, reverse

export ModPoly, poly, xpoly, field, coefficients, coeff, degree, deg
export LeadCoeff, ConstTerm, leading, setcoeff!, monic, ismonic
export DivRem, GCD, XGCD, derivative, diff, evaluate, compose
export LeftShift, RightShift, truncate, MulTrunc, InvTrunc
export MulMod, InvMod, InvModStatus, PowerMod, CompMod
export random_poly, BuildFromRoots, interpolate
export SquareFreeDecomp, DistinctDegreeDecomp, EqualDegreeSplit
export PolyFactorization, factor, reconstruct, is_irreducible, BuildIrred
export FindRoots, FindRoot

"""Dense univariate polynomial with coefficients stored in increasing degree order."""
struct ModPoly{T,F}
    _field::F
    _coeffs::Vector{T}
end

field(f::ModPoly) = f._field
@inline function field(f::ModPoly{PrimeField{P},DataType}) where {P}
    F = PrimeField{P}
    f._field === F || throw(ArgumentError("coefficient parent mismatch"))
    F
end
coefficients(f::ModPoly) = copy(f._coeffs)

function _trim!(c::Vector)
    while !isempty(c) && iszero(c[end])
        pop!(c)
    end
    c
end

function _make_poly(F::Type{PrimeField{P}}, c::AbstractVector) where {P}
    cc = Vector{PrimeField{P}}(c)
    ModPoly{PrimeField{P},DataType}(F, _trim!(cc))
end

function _make_poly(F, c::AbstractVector)
    # Empty comprehensions can otherwise degrade to Any[], especially when the
    # coefficient parent is itself a runtime extension-field object.
    T = typeof(zero(F))
    cc = Vector{T}(c)
    _trim!(cc)
    ModPoly{T,typeof(F)}(F, cc)
end

"""Construct a polynomial from low-to-high coefficients over parent `F`."""
function poly(F, input::AbstractVector)
    if isempty(input)
        T = typeof(zero(F))
        return _make_poly(F, Vector{T}())
    end
    c = [F(x) for x in input]
    _make_poly(F, c)
end
poly(F, c::Number) = poly(F, [c])
poly(F, f::ModPoly) = (field(f) == F || throw(ArgumentError("coefficient parent mismatch")); f)

"""Return the indeterminate `x` over `F`."""
xpoly(F) = poly(F, [zero(F), one(F)])

zero(f::ModPoly) = _make_poly(field(f), eltype(f._coeffs)[])
one(f::ModPoly) = poly(field(f), one(field(f)))
iszero(f::ModPoly) = isempty(f._coeffs)
isone(f::ModPoly) = degree(f) == 0 && isone(LeadCoeff(f))

"""Polynomial degree; the zero polynomial has degree `-1`."""
degree(f::ModPoly) = length(f._coeffs) - 1
deg(f::ModPoly) = degree(f)

"""Return coefficient `i`; out-of-range coefficients are zero."""
function coeff(f::ModPoly, i::Integer)
    i >= 0 || throw(DomainError(i, "a coefficient index must be nonnegative"))
    i < length(f._coeffs) ? f._coeffs[Int(i) + 1] : zero(field(f))
end
LeadCoeff(f::ModPoly) = iszero(f) ? zero(field(f)) : f._coeffs[end]
ConstTerm(f::ModPoly) = coeff(f, 0)
leading(f::ModPoly) = LeadCoeff(f)

"""Set a coefficient in place and normalize trailing zeros."""
function setcoeff!(f::ModPoly, i::Integer, value)
    i >= 0 || throw(DomainError(i, "a coefficient index must be nonnegative"))
    j = Int(i) + 1
    while length(f._coeffs) < j
        push!(f._coeffs, zero(field(f)))
    end
    f._coeffs[j] = field(f)(value)
    _trim!(f._coeffs)
    f
end

==(f::ModPoly, g::ModPoly) = field(f) == field(g) && f._coeffs == g._coeffs
isequal(f::ModPoly, g::ModPoly) = f == g

function show(io::IO, f::ModPoly)
    print(io, "ModPoly(", field(f), ", [")
    for (i, c) in enumerate(f._coeffs)
        i > 1 && print(io, ", ")
        show(io, c)
    end
    print(io, "])" )
end

function _check_poly_parents(f::ModPoly, g::ModPoly)
    field(f) == field(g) || throw(ArgumentError("polynomial coefficient parents differ"))
    field(f)
end

+(f::ModPoly, g::ModPoly) = begin
    F = _check_poly_parents(f, g)
    n = max(length(f._coeffs), length(g._coeffs))
    c = [zero(F) for _ in 1:n]
    for i in 1:n
        c[i] = coeff(f, i - 1) + coeff(g, i - 1)
    end
    _make_poly(F, c)
end

-(f::ModPoly, g::ModPoly) = begin
    F = _check_poly_parents(f, g)
    n = max(length(f._coeffs), length(g._coeffs))
    c = [coeff(f, i - 1) - coeff(g, i - 1) for i in 1:n]
    _make_poly(F, c)
end
-(f::ModPoly) = _make_poly(field(f), [-c for c in f._coeffs])

function poly_mul_schoolbook(f::ModPoly, g::ModPoly)
    F = _check_poly_parents(f, g)
    (iszero(f) || iszero(g)) && return zero(f)
    c = fill(zero(F), degree(f) + degree(g) + 1)
    for i in 0:degree(f), j in 0:degree(g)
        c[i + j + 1] += f._coeffs[i + 1] * g._coeffs[j + 1]
    end
    _make_poly(F, c)
end

+(f::ModPoly, c::Number) = f + poly(field(f), c)
+(c::Number, f::ModPoly) = poly(field(f), c) + f
-(f::ModPoly, c::Number) = f - poly(field(f), c)
-(c::Number, f::ModPoly) = poly(field(f), c) - f
*(f::ModPoly, c::Number) = _make_poly(field(f), [a * field(f)(c) for a in f._coeffs])
*(c::Number, f::ModPoly) = f * c

function ^(f::ModPoly, e::Integer)
    e >= 0 || throw(DomainError(e, "ordinary polynomial powers must be nonnegative"))
    k = BigInt(e)
    result = one(f)
    base = f
    while k > 0
        isodd(k) && (result *= base)
        k >>= 1
        k > 0 && (base *= base)
    end
    result
end

"""Scale a polynomial by `x^n`."""
LeftShift(f::ModPoly, n::Integer) = begin
    n >= 0 || throw(DomainError(n, "the shift must be nonnegative"))
    iszero(f) ? zero(f) : _make_poly(field(f), vcat([zero(field(f)) for _ in 1:Int(n)], f._coeffs))
end

"""Discard coefficients below degree `n`."""
RightShift(f::ModPoly, n::Integer) = begin
    n >= 0 || throw(DomainError(n, "the shift must be nonnegative"))
    n > degree(f) ? zero(f) : _make_poly(field(f), f._coeffs[(Int(n) + 1):end])
end
<<(f::ModPoly, n::Integer) = LeftShift(f, n)
>>(f::ModPoly, n::Integer) = RightShift(f, n)

"""Polynomial long division. Return `(quotient, remainder)`."""
function _divrem_buffers(f::ModPoly, g::ModPoly, want_quotient::Bool)
    F = _check_poly_parents(f, g)
    iszero(g) && throw(DivideError())
    df, dg = degree(f), degree(g)
    T = typeof(zero(F))
    df < dg && return (T[], copy(f._coeffs))

    q = want_quotient ? fill(zero(F), df - dg + 1) : T[]
    r = copy(f._coeffs)
    lead = LeadCoeff(g)
    lead_inv = isone(lead) ? lead : inv(lead)
    @inbounds while length(r) > dg
        k = length(r) - dg - 1
        t = r[end] * lead_inv
        want_quotient && (q[k + 1] = t)
        # Update the remainder in its existing coefficient buffer. No shifted
        # polynomial, scalar multiple, or subtraction temporary is allocated.
        for j in 0:(dg - 1)
            r[k + j + 1] -= t * g._coeffs[j + 1]
        end
        pop!(r)
        _trim!(r)
    end
    q, r
end

function DivRem(f::ModPoly, g::ModPoly)
    F = _check_poly_parents(f, g)
    q, r = _divrem_buffers(f, g, true)
    (_make_poly(F, q), _make_poly(F, r))
end

function div(f::ModPoly, g::ModPoly)
    q, _ = _divrem_buffers(f, g, true)
    _make_poly(field(f), q)
end

function rem(f::ModPoly, g::ModPoly)
    _, r = _divrem_buffers(f, g, false)
    _make_poly(field(f), r)
end

/(f::ModPoly, g::ModPoly) = div(f, g)
function /(f::ModPoly, c::Number)
    F = field(f)
    iszero(f) && return zero(f)
    scalar = F(c)
    scale = isone(scalar) ? scalar : inv(scalar)
    _make_poly(F, [a * scale for a in f._coeffs])
end

"""Divide by the leading coefficient; zero is returned unchanged."""
monic(f::ModPoly) = iszero(f) ? f : f / LeadCoeff(f)
ismonic(f::ModPoly) = !iszero(f) && isone(LeadCoeff(f))

"""Monic Euclidean greatest common divisor."""
function gcd(f::ModPoly, g::ModPoly)
    _check_poly_parents(f, g)
    a, b = f, g
    while !iszero(b)
        a, b = b, rem(a, b)
    end
    monic(a)
end
GCD(f::ModPoly, g::ModPoly) = gcd(f, g)

"""Extended Euclidean algorithm returning `(d,s,t)` with `d=s*f+t*g`."""
function XGCD(f::ModPoly, g::ModPoly)
    _check_poly_parents(f, g)
    old_r, r = f, g
    old_s, s = one(f), zero(f)
    old_t, t = zero(f), one(f)
    while !iszero(r)
        q = div(old_r, r)
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    end
    if iszero(old_r)
        return old_r, old_s, old_t
    end
    scale = LeadCoeff(old_r)
    old_r / scale, old_s / scale, old_t / scale
end

"""Formal derivative."""
function derivative(f::ModPoly)
    F = field(f)
    degree(f) <= 0 && return zero(f)
    _make_poly(F, [F(i) * f._coeffs[i + 1] for i in 1:degree(f)])
end
diff(f::ModPoly) = derivative(f)

"""Horner evaluation at a scalar or a compatible field element."""
function evaluate(f::ModPoly, x)
    F = field(f)
    a = F(x)
    result = zero(F)
    for i in length(f._coeffs):-1:1
        result = result * a + f._coeffs[i]
    end
    result
end
evaluate(f::ModPoly, xs::AbstractVector) = [evaluate(f, x) for x in xs]

"""Polynomial composition `f(g(x))` by Horner's rule."""
function compose(f::ModPoly, g::ModPoly)
    _check_poly_parents(f, g)
    result = zero(f)
    for i in length(f._coeffs):-1:1
        result = result * g + f._coeffs[i]
    end
    result
end
CompMod(g::ModPoly, h::ModPoly, f::ModPoly) = rem(compose(g, h), f)

"""Keep only coefficients of degree less than `n`."""
function truncate(f::ModPoly, n::Integer)
    n >= 0 || throw(DomainError(n, "the truncation length must be nonnegative"))
    _make_poly(field(f), f._coeffs[1:min(Int(n), length(f._coeffs))])
end
MulTrunc(f::ModPoly, g::ModPoly, n::Integer) = truncate(f * g, n)

"""Newton iteration for the first `n` coefficients of `1/f`."""
function InvTrunc(f::ModPoly, n::Integer)
    n >= 0 || throw(DomainError(n, "the truncation length must be nonnegative"))
    n == 0 && return zero(f)
    iszero(ConstTerm(f)) && throw(DomainError(f, "the constant term must be invertible"))
    g = poly(field(f), inv(ConstTerm(f)))
    m = 1
    while m < n
        m2 = min(2m, Int(n))
        g = truncate(g * (poly(field(f), 2) - truncate(f, m2) * g), m2)
        m = m2
    end
    g
end

"""Multiplication followed by reduction modulo polynomial `m`."""
MulMod(f::ModPoly, g::ModPoly, m::ModPoly) = rem(f * g, m)

"""Inverse in the quotient ring modulo `m`."""
function InvMod(f::ModPoly, m::ModPoly)
    d, s, _ = XGCD(f, m)
    isone(d) || throw(DomainError(f, "polynomial inverse undefined"))
    rem(s, m)
end
function InvModStatus(f::ModPoly, m::ModPoly)
    d, s, _ = XGCD(f, m)
    isone(d) ? (true, rem(s, m)) : (false, zero(f))
end

"""Binary polynomial powering modulo `m`; negative exponents invert `f`."""
function PowerMod(f::ModPoly, e::Integer, m::ModPoly)
    base = rem(f, m)
    if e < 0
        base = InvMod(base, m)
    end
    magnitude = e isa Int ? unsigned(abs(e)) : abs(BigInt(e))
    _polynomial_power_nonnegative(base, magnitude, m)
end

_polynomial_power_nonnegative(base::ModPoly, k::Integer, m::ModPoly) =
    _polynomial_power_binary(base, k, m)

function _polynomial_power_binary(base::ModPoly, k::Integer, m::ModPoly)
    result = rem(one(base), m)
    while k > 0
        isodd(k) && (result = MulMod(result, base, m))
        k >>= 1
        k > 0 && (base = MulMod(base, base, m))
    end
    result
end

# For small quotients, use a fixed-size product buffer and recycle the
# reduced vectors across exponentiation steps.
function _small_polynomial_product!(out, work, a, b, normalized_modulus)
    fill!(work, zero(eltype(work)))
    d = length(out)
    @inbounds for i in 1:d
        for j in 1:d
            index = i + j - 1
            work[index] = muladd(a[j], b[i], work[index])
        end
    end
    @inbounds for shift in (d - 1):-1:1
        scale = -work[d + shift]
        for j in 1:d
            work[j + shift - 1] = muladd(scale, normalized_modulus[j], work[j + shift - 1])
        end
    end
    copyto!(out, 1, work, 1, d)
end

function _polynomial_power_nonnegative(base::ModPoly{T,DataType}, k::Integer,
                                       m::ModPoly{T,DataType}) where {T<:PrimeField}
    d = degree(m)
    if !(1 <= d <= 16 && k >= 4) || iszero(base)
        return _polynomial_power_binary(base, k, m)
    end
    lead = LeadCoeff(m)
    scale = isone(lead) ? lead : inv(lead)
    normalized_modulus = [m._coeffs[i] * scale for i in 1:d]
    value = fill(zero(T), d)
    copyto!(value, 1, base._coeffs, 1, length(base._coeffs))
    result = fill(zero(T), d)
    result[1] = one(T)
    temp = similar(result)
    work = Vector{T}(undef, 2d - 1)
    while !iszero(k)
        if isodd(k)
            _small_polynomial_product!(temp, work, result, value, normalized_modulus)
            result, temp = temp, result
        end
        k >>= 1
        if !iszero(k)
            _small_polynomial_product!(temp, work, value, value, normalized_modulus)
            value, temp = temp, value
        end
    end
    _make_poly(T, result)
end

"""Return a polynomial with `n` random coefficients, hence degree at most `n-1`."""
function random_poly(F, n::Integer; rng::AbstractRNG = Random.default_rng())
    n >= 0 || throw(DomainError(n, "the coefficient count must be nonnegative"))
    n == 0 && return poly(F, typeof(zero(F))[])
    _make_poly(F, [rand(rng, F) for _ in 1:Int(n)])
end

"""Build the monic polynomial whose roots are the supplied elements."""
function BuildFromRoots(F, roots::AbstractVector)
    f = poly(F, 1)
    x = xpoly(F)
    for r in roots
        f *= x - F(r)
    end
    f
end

"""Lagrange interpolation through distinct abscissae."""
function interpolate(F, xs::AbstractVector, ys::AbstractVector)
    length(xs) == length(ys) || throw(DimensionMismatch("x and y lengths differ"))
    isempty(xs) && throw(ArgumentError("interpolation needs at least one point"))
    length(unique(xs)) == length(xs) || throw(ArgumentError("interpolation abscissae must be distinct"))
    result = poly(F, 0)
    x = xpoly(F)
    for i in eachindex(xs)
        basis = poly(F, 1)
        denominator = one(F)
        for j in eachindex(xs)
            i == j && continue
            basis *= x - F(xs[j])
            denominator *= F(xs[i]) - F(xs[j])
        end
        result += basis * (F(ys[i]) / denominator)
    end
    result
end

# Factorization ---------------------------------------------------------------

_field_order(::Type{PrimeField{P}}) where {P} = BigInt(P)
_field_order(F::BigModRing) = modulus(F)
_field_characteristic(F) = BigInt(characteristic(F))
_validate_factor_parent(F) = nothing
function _validate_factor_parent(F::BigModRing)
    ProbPrime(modulus(F)) || throw(ArgumentError("polynomial factorization requires a prime-field parent"))
end

function _squarefree(f::ModPoly)
    iszero(f) && throw(DomainError(f, "zero has no square-free decomposition"))
    F = field(f)
    p = _field_characteristic(F)
    p > 1 || throw(ArgumentError("square-free decomposition needs positive characteristic"))
    f = monic(f)
    out = Tuple{typeof(f),Int}[]
    c = gcd(f, derivative(f))
    w = div(f, c)
    i = 1
    while !isone(w)
        y = gcd(w, c)
        fac = div(w, y)
        !isone(fac) && push!(out, (fac, i))
        w = y
        c = div(c, y)
        i += 1
    end
    if !isone(c)
        root = _pth_root_poly(c, p)
        for (fac, multiplicity) in _squarefree(root)
            push!(out, (fac, multiplicity * Int(p)))
        end
    end
    out
end

function _pth_root_poly(f::ModPoly, p::BigInt)
    F = field(f)
    c = [zero(F) for _ in 1:(Int(div(BigInt(degree(f)), p)) + 1)]
    for i in 0:degree(f)
        iszero(f._coeffs[i + 1]) && continue
        e = BigInt(i)
        e % p == 0 || throw(ArgumentError("polynomial is not a p-th power"))
        c[Int(div(e, p)) + 1] = f._coeffs[i + 1]
    end
    _make_poly(F, c)
end

"""Square-free factorization returning `(factor, multiplicity)` pairs."""
SquareFreeDecomp(f::ModPoly) = _squarefree(f)

"""Distinct-degree factorization of a square-free polynomial."""
function DistinctDegreeDecomp(f::ModPoly; rng::AbstractRNG = Random.default_rng())
    F = field(f)
    _validate_factor_parent(F)
    q = F isa Type && F <: PrimeField ? modulus(F) : _field_order(F)
    f = monic(f)
    degree(f) <= 1 && return degree(f) == 1 ? [(f, 1)] : Tuple{typeof(f),Int}[]
    x = xpoly(F)
    h = x
    rest = f
    out = Tuple{typeof(f),Int}[]
    d = 1
    while degree(rest) >= 2d
        h = PowerMod(h, q, rest)
        g = gcd(rest, h - x)
        if !isone(g)
            push!(out, (g, d))
            rest = div(rest, g)
            h = rem(h, rest)
        end
        d += 1
    end
    !isone(rest) && push!(out, (rest, degree(rest)))
    out
end

function _trace_power_two(f::ModPoly, d::Int, m::ModPoly)
    result = zero(f)
    term = f
    for _ in 1:d
        result = rem(result + term, m)
        term = PowerMod(term, 2, m)
    end
    result
end

"""Cantor–Zassenhaus equal-degree splitting."""
function EqualDegreeSplit(f::ModPoly, d::Integer; rng::AbstractRNG = Random.default_rng())
    d >= 1 || throw(DomainError(d, "the common irreducible degree must be positive"))
    F = field(f)
    _validate_factor_parent(F)
    q = _field_order(F)
    f = monic(f)
    degree(f) == d && return [f]
    degree(f) % d == 0 || throw(ArgumentError("the input is not an equal-degree product"))
    while true
        a = random_poly(F, degree(f); rng = rng)
        iszero(a) && continue
        g = gcd(f, a)
        if isone(g)
            if q == 2
                b = _trace_power_two(a, Int(d), f)
            else
                exponent = div(q^Int(d) - 1, 2)
                b = PowerMod(a, exponent, f) - one(f)
            end
            g = gcd(f, b)
        end
        if !isone(g) && g != f
            left = EqualDegreeSplit(g, d; rng = rng)
            right = EqualDegreeSplit(div(f, g), d; rng = rng)
            return vcat(left, right)
        end
    end
end

"""The unit and `(irreducible factor, multiplicity)` pairs of a polynomial."""
struct PolyFactorization{T,F}
    unit::T
    factors::Vector{Tuple{ModPoly{T,F},Int}}
end

function show(io::IO, fac::PolyFactorization)
    print(io, "PolyFactorization(", fac.unit, ", ")
    show(io, fac.factors)
    print(io, ")")
end

"""Reconstruct the original polynomial from its factorization."""
function reconstruct(fac::PolyFactorization)
    result = poly(parent(fac.unit), fac.unit)
    for (g, multiplicity) in fac.factors
        result *= g^multiplicity
    end
    result
end

"""Factor a nonzero polynomial over a finite field."""
function factor(f::ModPoly; rng::AbstractRNG = Random.default_rng())
    _factor_generic(f, rng)
end

function _factor_generic(f::ModPoly, rng::AbstractRNG)
    iszero(f) && throw(DomainError(f, "zero has no irreducible factorization"))
    F = field(f)
    _validate_factor_parent(F)
    _field_order(F) # Validate that a finite-field order is available.
    unit = LeadCoeff(f)
    factors = Tuple{typeof(f),Int}[]
    for (sf, multiplicity) in _squarefree(f)
        for (ddf, d) in DistinctDegreeDecomp(sf; rng = rng)
            for irreducible in EqualDegreeSplit(ddf, d; rng = rng)
                push!(factors, (monic(irreducible), multiplicity))
            end
        end
    end
    sort!(factors, by = pair -> (degree(pair[1]), pair[1]._coeffs == f._coeffs ? 1 : 0, pair[2]))
    PolyFactorization(unit, factors)
end

# Factor a monic root-free quartic into monic quadratics if possible.
# For (x²+a*x+b)(x²+c*x+d), c=A-a, b+d=B-a*c, and C=a*(b+d)+(c-a)*b.
function _quartic_pair(f::ModPoly{PrimeField{P},DataType}) where {P}
    F = field(f)
    D, C, B, A = f._coeffs
    for i in 0:P-1
        a = F(i)
        c = A - a
        total = B - a*c
        if a == c
            C == a*total || continue
            for j in 0:P-1
                b = F(j)
                d = total - b
                if b*d == D
                    return (poly(F, [b,a,one(F)]), poly(F, [d,c,one(F)]))
                end
            end
        else
            b = (C - a*total)/(c - a)
            d = total - b
            b*d == D && return (poly(F, [b,a,one(F)]), poly(F, [d,c,one(F)]))
        end
    end
    nothing
end

# On a tiny prime field, testing every field element avoids modular powers
# for the linear part. Synthetic division uses a private coefficient buffer.
function factor(f::ModPoly{PrimeField{P},DataType}; rng::AbstractRNG = Random.default_rng()) where {P}
    (P <= 127 && 1 <= degree(f) <= 16) || return _factor_generic(f, rng)
    F = field(f)
    unit = LeadCoeff(f)
    c = coefficients(monic(f))
    factors = Tuple{typeof(f),Int}[]
    for r in 0:P-1
        length(c) == 1 && break
        root = F(r)
        count = 0
        while length(c) > 1
            value = c[end]
            for i in length(c)-1:-1:1
                value = muladd(value, root, c[i])
            end
            iszero(value) || break
            carry = pop!(c)
            for i in length(c):-1:1
                old = c[i]
                c[i] = carry
                carry = muladd(carry, root, old)
            end
            count += 1
        end
        count > 0 && push!(factors, (poly(F, [-root, one(F)]), count))
    end
    if length(c) > 1
        rest = _make_poly(F, c)
        if degree(rest) <= 3
            push!(factors, (rest, 1)) # No root over F, hence irreducible.
        elseif degree(rest) == 4
            pair = _quartic_pair(rest)
            if isnothing(pair)
                push!(factors, (rest, 1))
            elseif pair[1] == pair[2]
                push!(factors, (pair[1], 2))
            else
                push!(factors, (pair[1], 1), (pair[2], 1))
            end
        else
            append!(factors, _factor_generic(rest, rng).factors)
        end
    end
    sort!(factors, by = pair -> (degree(pair[1]), pair[2]))
    PolyFactorization(unit, factors)
end

"""Rabin irreducibility test over a finite field."""
function is_irreducible(f::ModPoly)
    n = degree(f)
    n > 0 || return false
    F = field(f)
    _validate_factor_parent(F)
    q = _field_order(F)
    f = monic(f)
    n == 1 && return true
    x = xpoly(F)
    for (r, _) in FactorInt(n)
        exponent = q^div(n, r)
        h = PowerMod(x, exponent, f)
        isone(gcd(f, h - x)) || return false
    end
    PowerMod(x, q^n, f) - x == zero(f)
end

"""Generate a random monic irreducible polynomial of degree `n`."""
function BuildIrred(F, n::Integer; rng::AbstractRNG = Random.default_rng())
    n > 0 || throw(DomainError(n, "the degree must be positive"))
    while true
        f = random_poly(F, n; rng = rng) + LeftShift(poly(F, 1), n)
        is_irreducible(f) && return f
    end
end

"""Return all roots with multiplicity using polynomial factorization."""
function FindRoots(f::ModPoly; rng::AbstractRNG = Random.default_rng())
    fac = factor(f; rng = rng)
    F = field(f)
    out = eltype(f._coeffs)[]
    for (g, multiplicity) in fac.factors
        if degree(g) == 1
            root = -ConstTerm(g) / LeadCoeff(g)
            for _ in 1:multiplicity
                push!(out, root)
            end
        end
    end
    out
end

"""Return one root, or `nothing` when no root exists."""
function FindRoot(f::ModPoly; rng::AbstractRNG = Random.default_rng())
    roots = FindRoots(f; rng = rng)
    isempty(roots) ? nothing : roots[1]
end
