# Fast LLL with exact post-verification ----------------------------------------

_dot_float(a, b) = sum(a[i] * b[i] for i in eachindex(a, b))

function _gso_float(B::AbstractMatrix{<:Integer})
    n = size(B, 1)
    mu = zeros(Float64, n, n)
    norms = zeros(Float64, n)
    bstar = Vector{Vector{Float64}}()
    for i in 1:n
        row = Float64.(@view B[i, :])
        v = copy(row)
        mu[i, i] = 1.0
        for j in 1:(i - 1)
            iszero(norms[j]) && throw(ArgumentError("the input rows must be linearly independent"))
            mu[i, j] = _dot_float(row, bstar[j]) / norms[j]
            v .-= mu[i, j] .* bstar[j]
        end
        push!(bstar, v)
        norms[i] = _dot_float(v, v)
        iszero(norms[i]) && throw(ArgumentError("the input rows must be linearly independent"))
    end
    mu, norms
end

_nearest_float(x::Float64) = floor(Int, x + 0.5)

# An adjacent swap changes two orthogonal vectors. Update later rows'
# coordinates in that plane and leave the remaining coordinates in place.
function _lll_swap_gso!(mu, norms, k)
    j = k - 1
    a = mu[k, j]
    old_j, old_k = norms[j], norms[k]
    new_j = old_k + a^2 * old_j
    new_k = old_k * (old_j / new_j)
    isfinite(new_j) && isfinite(new_k) && new_j > 0 && new_k > 0 || return false
    new_a = a * (old_j / new_j)
    isfinite(new_a) || return false
    norms[j], norms[k] = new_j, new_k
    mu[k, j] = new_a
    @inbounds for i in 1:(j - 1)
        mu[j, i], mu[k, i] = mu[k, i], mu[j, i]
    end
    @inbounds for i in (k + 1):size(mu, 1)
        previous = mu[i, k]
        mu[i, k] = mu[i, j] - a * previous
        mu[i, j] = previous + new_a * mu[i, k]
    end
    true
end

function _lll_row_subtract!(A::Matrix{T}, k::Int, j::Int, q::Int) where {T<:Integer}
    @inbounds for i in axes(A, 2)
        A[k, i] -= q * A[j, i]
    end
end

function _lll_row_subtract!(A::Matrix{Int}, k::Int, j::Int, q::Int)
    @inbounds for i in axes(A, 2)
        A[k, i] = Base.Checked.checked_sub(A[k, i], Base.Checked.checked_mul(q, A[j, i]))
    end
end

function _lll_fast_reduce(A::AbstractMatrix{<:Integer}, d::Float64, maxiter::Int,
                          ::Type{T}) where {T<:Integer}
    n = size(A, 1)
    B = Matrix{T}(A)
    U = zeros(T, n, n)
    for i in 1:n
        U[i, i] = 1
    end
    n == 1 && return B, U
    mu, norms = _gso_float(B)
    k = 2
    iterations = 0
    swaps_since_refresh = 0
    while k <= n
        iterations += 1
        iterations <= maxiter || error("LLL iteration limit exceeded")
        for j in (k - 1):-1:1
            q = _nearest_float(mu[k, j])
            if q != 0
                _lll_row_subtract!(B, k, j, q)
                _lll_row_subtract!(U, k, j, q)
                for l in 1:j
                    mu[k, l] -= q * mu[j, l]
                end
            end
        end
        if norms[k] >= (d - mu[k, k - 1]^2) * norms[k - 1]
            k += 1
        else
            @inbounds for i in axes(B, 2)
                B[k - 1, i], B[k, i] = B[k, i], B[k - 1, i]
            end
            @inbounds for i in axes(U, 2)
                U[k - 1, i], U[k, i] = U[k, i], U[k - 1, i]
            end
            swaps_since_refresh += 1
            if swaps_since_refresh >= n || !_lll_swap_gso!(mu, norms, k)
                mu, norms = _gso_float(B)
                swaps_since_refresh = 0
            end
            k = max(k - 1, 2)
        end
    end
    B, U
end

function _lll_fast_core(A::AbstractMatrix{<:Integer}, d::Float64, maxiter::Int)
    if all(x -> typemin(Int) <= x <= typemax(Int), A)
        try
            return _lll_fast_reduce(A, d, maxiter, Int)
        catch err
            err isa OverflowError || rethrow()
        end
    end
    _lll_fast_reduce(A, d, maxiter, BigInt)
end

function _lll_transform_matches(A, B, U)
    Z = Matrix{BigInt}(A)
    acc, scratch, zero_value = BigInt(0), BigInt(0), BigInt(0)
    @inbounds for j in axes(B, 2), i in axes(B, 1)
        Base.GMP.MPZ.set!(acc, zero_value)
        for k in axes(Z, 1)
            Base.GMP.MPZ.mul!(scratch, U[i, k], Z[k, j])
            Base.GMP.MPZ.add!(acc, acc, scratch)
        end
        acc == B[i, j] || return false
    end
    true
end

function _lll_verify(A, B, U, d)
    try
        _lll_transform_matches(A, B, U) || return false
        abs(det_exact(U)) == 1 || return false
        gram = _gram_bareiss(B)
        for i in 1:size(B, 1), j in 1:(i - 1)
            abs(2 * gram[i, j]) <= gram[j, j] || return false
        end
        num, den = numerator(d), denominator(d)
        prior = BigInt(1)
        for k in 2:size(B, 1)
            pivot = gram[k - 1, k - 1]
            # Cross-multiply the Lovász inequality by the positive Gram pivots.
            den * (gram[k, k] * prior + gram[k, k - 1]^2) >= num * pivot^2 || return false
            prior = pivot
        end
        true
    catch
        false
    end
end

# Exact determinant up to sign. Row exchanges do not affect the predicate.
function _lll_word_unimodular(U::Matrix{Int})
    M = copy(U)
    n = size(M, 1)
    previous = 1
    @inbounds for k in 1:n-1
        if iszero(M[k,k])
            row = k + 1
            while row <= n && iszero(M[row,k])
                row += 1
            end
            row > n && return false
            for j in k:n
                M[k,j], M[row,j] = M[row,j], M[k,j]
            end
        end
        pivot = M[k,k]
        for j in k+1:n, i in k+1:n
            left = Base.Checked.checked_mul(M[i,j], pivot)
            right = Base.Checked.checked_mul(M[i,k], M[k,j])
            M[i,j] = div(Base.Checked.checked_sub(left, right), previous)
        end
        previous = pivot
    end
    M[n,n] == 1 || M[n,n] == -1
end

# The caller bounds each input by 2^20 and both dimensions by 128:
# dot products fit signed Int64 (at most 2^47 in magnitude).
function _lll_word_verify(A, B, U, d)
    n, m = size(B)
    @inbounds for j in 1:m, i in 1:n
        total = 0
        for k in 1:n
            total += U[i,k] * A[k,j]
        end
        total == B[i,j] || return false
    end
    _lll_word_unimodular(U) || return false
    gram = Matrix{Int128}(undef, n, n)
    @inbounds for j in 1:n, i in j:n
        total = 0
        for k in 1:m
            total += B[i,k] * B[j,k]
        end
        gram[i,j] = total
    end
    mul = Base.Checked.checked_mul
    sub = Base.Checked.checked_sub
    add = Base.Checked.checked_add
    num, den = Int128(numerator(d)), Int128(denominator(d))
    previous = Int128(1)
    @inbounds for j in 1:n
        pivot = gram[j,j]
        pivot > 0 || return false
        for i in j+1:n
            # Equivalent to 2*abs(mu_numerator) <= pivot, without overflow.
            half = pivot >> 1
            -half <= gram[i,j] <= half || return false
        end
        for k in j+1:n, i in k:n
            gram[i,k] = div(sub(mul(gram[i,k], pivot), mul(gram[i,j], gram[k,j])), previous)
        end
        if j < n
            off = gram[j+1,j]
            left = mul(den, add(mul(gram[j+1,j+1], previous), mul(off,off)))
            left >= mul(num, mul(pivot,pivot)) || return false
        end
        previous = pivot
    end
    true
end

function _lll_verify_reduction(A, B::Matrix{Int}, U::Matrix{Int}, d, bigB, bigU)
    small(x) = -(1<<20) <= x <= (1<<20)
    if Sys.WORD_SIZE == 64 && 1 <= size(B,1) <= 12 && 1 <= size(B,2) <= 128 &&
       all(small, B) && all(small, U) &&
       numerator(d) <= typemax(Int128) && denominator(d) <= typemax(Int128)
        if all(small, A)
            # Check before converting A: an oversized original matrix must use
            # the generic BigInt verifier rather than fail during conversion.
            Z = Matrix{Int}(A)
            try
                return _lll_word_verify(Z, B, U, d)
            catch err
                err isa OverflowError || err isa DivideError || rethrow()
            end
        end
    end
    _lll_verify(A, bigB, bigU, d)
end
_lll_verify_reduction(A, B, U, d, bigB, bigU) = _lll_verify(A, bigB, bigU, d)

"""
    lll(A; delta=3//4, mode=:fast, verify=true, maxiter=100000)

LLL with exact integer row operations. `:fast` uses floating GSO only for
reduction decisions and verifies the returned basis with exact rational GSO.
`:exact` preserves the original exact-rational implementation.
"""
function lll(A::AbstractMatrix{<:Integer}; delta = BigInt(3) // BigInt(4),
             mode::Symbol = :fast, verify::Bool = true, maxiter::Integer = 100000)
    d = _lll_delta(delta)
    maxiter > 0 || throw(DomainError(maxiter, "maxiter must be positive"))
    if mode == :exact
        B, U = _lll_exact_legacy(A; delta = d, maxiter = maxiter)
        return (basis = B, transform = U, mode = :exact, verified = true)
    end
    mode == :fast || throw(ArgumentError("unknown LLL mode: $mode"))
    size(A, 1) == 1 && return (basis = BigInt.(A), transform = ones(BigInt, 1, 1), mode = :fast, verified = true)
    try
        B, U = _lll_fast_core(A, Float64(d), Int(maxiter))
        bigB, bigU = convert(Matrix{BigInt}, B), convert(Matrix{BigInt}, U)
        if !verify || _lll_verify_reduction(A, B, U, d, bigB, bigU)
            return (basis = bigB, transform = bigU, mode = :fast, verified = verify)
        end
    catch err
        err isa InterruptException && rethrow()
    end
    B, U = _lll_exact_legacy(A; delta = d, maxiter = maxiter)
    (basis = B, transform = U, mode = :exact, verified = true)
end
