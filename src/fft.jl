# Generalized FFT convolution with exact coefficient reconstruction ----------

export fft_convolve, FFT_PRIMES

# These primes support large radix-two transforms. The implementation selects
# only primes compatible with the requested transform length.
const FFT_PRIMES = (
    998244353,   # 119 * 2^23 + 1
    469762049,   # 7 * 2^26 + 1
    985661441,   # 235 * 2^22 + 1
    943718401,   # 225 * 2^22 + 1
    1004535809,  # 479 * 2^21 + 1
    2013265921,  # 15 * 2^27 + 1
)

# Auxiliary NTT primes are runtime values. Use a function barrier so the
# selected modulus becomes a compile-time parameter inside the actual kernel.
function _fft_ntt_primefield(a, b, q::Int)
    _fft_ntt_primefield_typed(a, b, PrimeField(q))
end

function _fft_ntt_primefield_typed(a, b, ::Type{PrimeField{Q}}) where {Q}
    A = Vector{PrimeField{Q}}(undef, length(a))
    B = Vector{PrimeField{Q}}(undef, length(b))
    @inbounds for i in eachindex(a)
        A[i] = _raw(Val{Q}(), Int(mod(rep(a[i]), Q)))
    end
    @inbounds for i in eachindex(b)
        B[i] = _raw(Val{Q}(), Int(mod(rep(b[i]), Q)))
    end
    # Q is one of the validated FFT primes. Raw construction is therefore safe,
    # and avoids repeating primality checks for every coefficient.
    [rep(x) for x in ntt_convolve(A, B)]
end

function _fft_bound(a, b)
    isempty(a) && return BigInt(0)
    isempty(b) && return BigInt(0)
    ma = maximum(abs(rep(x)) for x in a)
    mb = maximum(abs(rep(x)) for x in b)
    BigInt(ma) * BigInt(mb) * min(length(a), length(b))
end

function _fft_direct_supported(F, total)
    F == GF2 && return false
    P = modulus(F)
    n = nextpow(2, total)
    (P - 1) % n == 0 || return false
    F isa BigModRing ? ProbPrime(P) : true
end

function _fft_prime_plan(bound::BigInt, total::Int)
    n = nextpow(2, total)
    candidates = [q for q in FFT_PRIMES if (q - 1) % n == 0]
    isempty(candidates) && throw(ArgumentError("no FFT prime supports the required transform length $n"))
    moduli = Int[]
    product = BigInt(1)
    for q in candidates
        push!(moduli, q)
        product *= q
        product > bound && break
    end
    product > bound || throw(ArgumentError("the coefficient bound exceeds the available FFT prime product"))
    moduli
end

function _crt_integers(residues::Vector{BigInt}, moduli::Vector{Int})
    x = BigInt(0)
    M = BigInt(1)
    for (r, q) in zip(residues, moduli)
        t = mod((r - x) * invmod(M, q), q)
        x = mod(x + M * t, M * q)
        M *= q
    end
    x, M
end

"""
    fft_convolve(a, b)

Compute an exact integer convolution using FFTs over one or more NTT primes and
reconstruct the result in the input scalar parent. The coefficient bound is
computed from the inputs, so reconstruction is exact whenever the selected
prime product exceeds that bound.
"""
function fft_convolve(a::AbstractVector, b::AbstractVector)
    isempty(a) && return copy(a)
    isempty(b) && return copy(b)
    length(a) == 0 && return copy(a)
    parent(a[1]) == parent(b[1]) || throw(ArgumentError("vector parents differ"))
    F = parent(a[1])
    total = length(a) + length(b) - 1
    _fft_direct_supported(F, total) && return ntt_convolve(a, b)
    bound = _fft_bound(a, b)
    moduli = _fft_prime_plan(bound, total)
    residues = Vector{Vector{Int}}()
    for q in moduli
        push!(residues, _fft_ntt_primefield(a, b, q))
    end
    F = parent(a[1])
    if length(moduli) == 1
        return _fft_reconstruct(F, residues[1])
    end
    F.(_crt_coefficients(residues, moduli, total))
end

function _crt_coefficients(residues::Vector{<:AbstractVector{<:Integer}}, moduli::Vector{Int}, total::Int)
    out = zeros(BigInt, total)
    M = BigInt(1)
    for (k, q) in enumerate(moduli)
        # Compute the inverse once per prime, not once per coefficient.
        factor = invmod(M, q)
        for i in 1:total
            t = mod((residues[k][i] - out[i]) * factor, q)
            out[i] = mod(out[i] + M * t, M * q)
        end
        M *= q
    end
    out
end

function _fft_reconstruct(::Type{PrimeField{P}}, residues::AbstractVector{<:Integer}) where {P}
    # A single auxiliary prime proves that the nonnegative integer coefficients
    # are already exact. Keep reconstruction in machine integers on this path.
    [_raw(Val{P}(), Int(mod(x, P))) for x in residues]
end

function _fft_reconstruct(F::BigModRing, residues::AbstractVector{<:Integer})
    F.(residues)
end

function _fft_reconstruct(::Type{GF2}, residues::AbstractVector{<:Integer})
    GF2.(residues)
end

function fft_convolve(a::Vector{PrimeField{P}}, b::Vector{PrimeField{P}}) where {P}
    isempty(a) && return copy(a)
    isempty(b) && return copy(b)
    total = length(a) + length(b) - 1
    _fft_direct_supported(PrimeField{P}, total) && return ntt_convolve(a, b)
    bound = _fft_bound(a, b)
    moduli = _fft_prime_plan(bound, total)
    residues = Vector{Vector{Int}}()
    for q in moduli
        push!(residues, _fft_ntt_primefield(a, b, q))
    end
    if length(moduli) == 1
        return _fft_reconstruct(PrimeField{P}, residues[1])
    end
    _fft_reconstruct(PrimeField{P}, _crt_coefficients(residues, moduli, total))
end

function fft_convolve(a::Vector{ModInt}, b::Vector{ModInt})
    isempty(a) && return copy(a)
    isempty(b) && return copy(b)
    F = parent(a[1])
    parent(b[1]) == F || throw(ArgumentError("vector parents differ"))
    total = length(a) + length(b) - 1
    _fft_direct_supported(F, total) && return ntt_convolve(a, b)
    bound = _fft_bound(a, b)
    moduli = _fft_prime_plan(bound, total)
    residues = Vector{Vector{Int}}()
    for q in moduli
        push!(residues, _fft_ntt_primefield(a, b, q))
    end
    if length(moduli) == 1
        return _fft_reconstruct(F, residues[1])
    end
    _fft_reconstruct(F, _crt_coefficients(residues, moduli, total))
end

function fft_convolve(a::Vector{GF2}, b::Vector{GF2})
    isempty(a) && return copy(a)
    isempty(b) && return copy(b)
    total = length(a) + length(b) - 1
    _fft_direct_supported(GF2, total) && return ntt_convolve(a, b)
    bound = _fft_bound(a, b)
    moduli = _fft_prime_plan(bound, total)
    residues = Vector{Vector{Int}}()
    for q in moduli
        push!(residues, _fft_ntt_primefield(a, b, q))
    end
    if length(moduli) == 1
        return _fft_reconstruct(GF2, residues[1])
    end
    _fft_reconstruct(GF2, _crt_coefficients(residues, moduli, total))
end
