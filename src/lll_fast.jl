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
            B, U = _lll_fast_reduce(A, d, maxiter, Int)
            return BigInt.(B), BigInt.(U)
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
        if !verify || _lll_verify(A, B, U, d)
            return (basis = B, transform = U, mode = :fast, verified = verify)
        end
    catch err
        err isa InterruptException && rethrow()
    end
    B, U = _lll_exact_legacy(A; delta = d, maxiter = maxiter)
    (basis = B, transform = U, mode = :exact, verified = true)
end
