# Optional polynomial multiplication acceleration ------------------------------

export set_polynomial_acceleration!, polynomial_acceleration, poly_mul

const POLY_ACCELERATION = Ref(true)
const KARATSUBA_THRESHOLD = 64
const NTT_THRESHOLD = 256

"""Enable or disable degree-based polynomial multiplication acceleration."""
function set_polynomial_acceleration!(enabled::Bool = true)
    POLY_ACCELERATION[] = enabled
    enabled
end

"""Return whether degree-based polynomial acceleration is enabled."""
polynomial_acceleration() = POLY_ACCELERATION[]

function _schoolbook_coeff(a, b, z)
    (isempty(a) || isempty(b)) && return typeof(z)[]
    c = fill(z, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i + j - 1] += a[i] * b[j]
    end
    c
end

function _karatsuba_coeff(a, b, z)
    n = max(length(a), length(b))
    n <= KARATSUBA_THRESHOLD && return _schoolbook_coeff(a, b, z)
    n = nextpow(2, n)
    a = vcat(a, fill(z, n - length(a)))
    b = vcat(b, fill(z, n - length(b)))
    m = n >> 1
    a0, a1 = a[1:m], a[(m + 1):n]
    b0, b1 = b[1:m], b[(m + 1):n]
    z0 = _karatsuba_coeff(a0, b0, z)
    z2 = _karatsuba_coeff(a1, b1, z)
    z1 = _karatsuba_coeff(a0 + a1, b0 + b1, z)
    for i in eachindex(z0)
        z1[i] -= z0[i]
        z1[i] -= z2[i]
    end
    out = fill(z, 2n - 1)
    for i in eachindex(z0)
        out[i] += z0[i]
        out[i + n] += z2[i]
        out[i + m] += z1[i]
    end
    out
end

function _poly_mul_karatsuba(f::ModPoly, g::ModPoly)
    F = _check_poly_parents(f, g)
    _make_poly(F, _karatsuba_coeff(copy(f._coeffs), copy(g._coeffs), zero(F)))
end

function _ntt_length_supported(F, total)
    F == GF2 && return false
    P = modulus(F)
    n = nextpow(2, total)
    (P - 1) % n == 0 || return false
    F isa ZZPField ? ProbPrime(P) : true
end

function _poly_mul_ntt(f::ModPoly, g::ModPoly)
    F = _check_poly_parents(f, g)
    total = degree(f) + degree(g) + 1
    _ntt_length_supported(F, total) || throw(ArgumentError("the selected modulus does not support the required NTT length"))
    _make_poly(F, ntt_convolve(f._coeffs, g._coeffs))
end

_fft_parent_supported(F) = (F isa Type && F <: zzp) || F isa ZZPField || F == GF2

function _poly_mul_fft(f::ModPoly, g::ModPoly)
    F = _check_poly_parents(f, g)
    _fft_parent_supported(F) || throw(ArgumentError("FFT convolution requires scalar modular coefficients"))
    _make_poly(F, fft_convolve(f._coeffs, g._coeffs))
end

function _poly_mul_auto(f::ModPoly, g::ModPoly)
    n = max(length(f._coeffs), length(g._coeffs))
    n <= KARATSUBA_THRESHOLD && return poly_mul_schoolbook(f, g)
    n < NTT_THRESHOLD && return _poly_mul_karatsuba(f, g)
    _fft_parent_supported(field(f)) || return _poly_mul_karatsuba(f, g)
    try
        return _poly_mul_fft(f, g)
    catch err
        err isa ArgumentError || rethrow()
        return _poly_mul_karatsuba(f, g)
    end
end

"""
    poly_mul(f, g; method=:default)

Multiply two dense polynomials. `method` can be `:default`, `:schoolbook`,
`:karatsuba`, `:ntt`, `:fft`, or `:auto`. The default preserves the original schoolbook
algorithm unless acceleration has been disabled explicitly. Acceleration is enabled by default.
"""
function poly_mul(f::ModPoly, g::ModPoly; method::Symbol = :default)
    method == :schoolbook && return poly_mul_schoolbook(f, g)
    method == :karatsuba && return _poly_mul_karatsuba(f, g)
    method == :ntt && return _poly_mul_ntt(f, g)
    method == :fft && return _poly_mul_fft(f, g)
    (method == :default || method == :auto) || throw(ArgumentError("unknown polynomial multiplication method: $method"))
    if method == :auto || POLY_ACCELERATION[]
        return _poly_mul_auto(f, g)
    end
    poly_mul_schoolbook(f, g)
end

*(f::ModPoly, g::ModPoly) = poly_mul(f, g)
