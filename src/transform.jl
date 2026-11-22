# Number-theoretic transforms --------------------------------------------------

export primitive_root, root_of_unity, ntt, intt, convolve, ntt_convolve

"""Return a primitive element of a supported prime field."""
primitive_root(F) = F(PrimitiveRoot(modulus(F)))

"""Return a primitive `n`-th root of unity in prime field `F`."""
function root_of_unity(n::Integer, F)
    n > 0 || throw(DomainError(n, "the root order must be positive"))
    P = modulus(F)
    (P - 1) % n == 0 || throw(ArgumentError("n must divide p-1"))
    primitive_root(F)^div(P - 1, n)
end

# Specialized in-place NTT for machine-word prime fields. The generic methods
# below remain available for other parents, while this path avoids dynamic
# field dispatch in every butterfly operation.
function _ntt_primefield_root!(A::Vector{PrimeField{P}}, omega::PrimeField{P}, inverse::Bool) where {P}
    n = length(A)
    n > 0 || throw(ArgumentError("NTT input must be nonempty"))
    n & (n - 1) == 0 || throw(ArgumentError("the NTT length must be a power of two"))
    n == 1 && return A
    F = PrimeField{P}

    j = 0
    @inbounds for i in 1:(n - 1)
        bit = n >> 1
        while (j & bit) != 0
            j = xor(j, bit)
            bit >>= 1
        end
        j = xor(j, bit)
        if i < j
            A[i + 1], A[j + 1] = A[j + 1], A[i + 1]
        end
    end

    one_element = _raw(Val{P}(), 1)
    len = 2
    while len <= n
        wlen = _raw(Val{P}(), PowerMod(rep(omega), n ÷ len, P))
        half = len >> 1
        @inbounds for start in 1:len:n
            w = one_element
            for offset in 0:(half - 1)
                i = start + offset
                j = i + half
                u = rep(A[i])
                v = _mulmod_raw(rep(A[j]), rep(w), P)
                A[i] = _raw(Val{P}(), _addmod_raw(u, v, P))
                A[j] = _raw(Val{P}(), _submod_raw(u, v, P))
                w = _raw(Val{P}(), _mulmod_raw(rep(w), rep(wlen), P))
            end
        end
        len <<= 1
    end
    if inverse
        scale = InvMod(n, P)
        @inbounds for i in eachindex(A)
            A[i] = _raw(Val{P}(), _mulmod_raw(rep(A[i]), scale, P))
        end
    end
    A
end

function _ntt_primefield!(A::Vector{PrimeField{P}}, inverse::Bool = false) where {P}
    length(A) == 1 && return A
    omega = inverse ? inv(root_of_unity(length(A), PrimeField{P})) :
                      root_of_unity(length(A), PrimeField{P})
    _ntt_primefield_root!(A, omega, inverse)
end

ntt(a::Vector{PrimeField{P}}) where {P} = _ntt_primefield!(collect(a), false)
intt(a::Vector{PrimeField{P}}) where {P} = _ntt_primefield!(collect(a), true)

function ntt_convolve(a::Vector{PrimeField{P}}, b::Vector{PrimeField{P}}) where {P}
    isempty(a) && return copy(a)
    isempty(b) && return copy(b)
    total = length(a) + length(b) - 1
    n = nextpow(2, total)
    F = PrimeField{P}
    z = _raw(Val{P}(), 0)
    A = fill(z, n)
    B = fill(z, n)
    copyto!(A, 1, a, 1, length(a))
    copyto!(B, 1, b, 1, length(b))
    # One primitive root computation is shared by both forward transforms and
    # the inverse transform of this convolution.
    omega = root_of_unity(n, F)
    _ntt_primefield_root!(A, omega, false)
    _ntt_primefield_root!(B, omega, false)
    @inbounds for i in 1:n
        A[i] = _raw(Val{P}(), _mulmod_raw(rep(A[i]), rep(B[i]), P))
    end
    _ntt_primefield_root!(A, inv(omega), true)
    A[1:total]
end

"""In-place iterative radix-two Cooley–Tukey NTT on a copied vector."""
function ntt(a::AbstractVector)
    n = length(a)
    n > 0 || throw(ArgumentError("NTT input must be nonempty"))
    n & (n - 1) == 0 || throw(ArgumentError("the NTT length must be a power of two"))
    A = collect(a)
    n == 1 && return A
    F = parent(A[1])
    omega = root_of_unity(n, F)

    j = 0
    for i in 1:(n - 1)
        bit = n >> 1
        while (j & bit) != 0
            j = xor(j, bit)
            bit >>= 1
        end
        j = xor(j, bit)
        if i < j
            A[i + 1], A[j + 1] = A[j + 1], A[i + 1]
        end
    end

    one_element = one(F)
    len = 2
    while len <= n
        wlen = omega^(n ÷ len)
        half = len >> 1
        for start in 1:len:n
            w = one_element
            for offset in 0:(half - 1)
                u = A[start + offset]
                v = A[start + offset + half] * w
                A[start + offset] = u + v
                A[start + offset + half] = u - v
                w *= wlen
            end
        end
        len <<= 1
    end
    A
end

"""Inverse NTT using the inverse root and division by the transform length."""
function intt(a::AbstractVector)
    n = length(a)
    n > 0 || throw(ArgumentError("NTT input must be nonempty"))
    n & (n - 1) == 0 || throw(ArgumentError("the NTT length must be a power of two"))
    n == 1 && return collect(a)
    F = parent(a[1])
    omega_inv = inv(root_of_unity(n, F))
    A = collect(a)

    j = 0
    for i in 1:(n - 1)
        bit = n >> 1
        while (j & bit) != 0
            j = xor(j, bit)
            bit >>= 1
        end
        j = xor(j, bit)
        if i < j
            A[i + 1], A[j + 1] = A[j + 1], A[i + 1]
        end
    end
    one_element = one(F)
    len = 2
    while len <= n
        wlen = omega_inv^(n ÷ len)
        half = len >> 1
        for start in 1:len:n
            w = one_element
            for offset in 0:(half - 1)
                u = A[start + offset]
                v = A[start + offset + half] * w
                A[start + offset] = u + v
                A[start + offset + half] = u - v
                w *= wlen
            end
        end
        len <<= 1
    end
    scale = inv(F(n))
    A .* scale
end

"""Schoolbook convolution for arbitrary supported field vectors."""
function convolve(a::AbstractVector, b::AbstractVector)
    isempty(a) && return copy(a)
    isempty(b) && return copy(b)
    F = parent(a[1])
    parent(b[1]) == F || throw(ArgumentError("vector parents differ"))
    c = [zero(F) for _ in 1:(length(a) + length(b) - 1)]
    for i in eachindex(a), j in eachindex(b)
        c[i + j - 1] += a[i] * b[j]
    end
    c
end

"""Zero-padded cyclic NTT convolution, trimmed to linear-convolution length."""
function ntt_convolve(a::AbstractVector, b::AbstractVector)
    isempty(a) && return copy(a)
    isempty(b) && return copy(b)
    F = parent(a[1])
    parent(b[1]) == F || throw(ArgumentError("vector parents differ"))
    total = length(a) + length(b) - 1
    n = nextpow(2, total)
    A = fill(zero(F), n)
    B = fill(zero(F), n)
    A[1:length(a)] = a
    B[1:length(b)] = b
    C = intt(ntt(A) .* ntt(B))
    C[1:total]
end
