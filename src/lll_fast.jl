# Fast LLL with exact post-verification ----------------------------------------

_dot_float(a, b) = sum(a[i] * b[i] for i in eachindex(a, b))

function _gso_float(B::AbstractMatrix{BigInt})
    n = size(B, 1)
    mu = zeros(Float64, n, n)
    norms = zeros(Float64, n)
    bstar = Vector{Vector{Float64}}()
    for i in 1:n
        v = Float64.(B[i, :])
        mu[i, i] = 1.0
        for j in 1:(i - 1)
            iszero(norms[j]) && throw(ArgumentError("the input rows must be linearly independent"))
            mu[i, j] = _dot_float(Float64.(B[i, :]), bstar[j]) / norms[j]
            v -= mu[i, j] .* bstar[j]
        end
        push!(bstar, v)
        norms[i] = _dot_float(v, v)
        iszero(norms[i]) && throw(ArgumentError("the input rows must be linearly independent"))
    end
    mu, norms
end

_nearest_float(x::Float64) = floor(Int, x + 0.5)

function _lll_fast_core(A::AbstractMatrix{<:Integer}, d::Float64, maxiter::Int)
    n = size(A, 1)
    B = BigInt.(A)
    U = zeros(BigInt, n, n)
    for i in 1:n
        U[i, i] = 1
    end
    n == 1 && return B, U
    mu, norms = _gso_float(B)
    k = 2
    iterations = 0
    while k <= n
        iterations += 1
        iterations <= maxiter || error("LLL iteration limit exceeded")
        for j in (k - 1):-1:1
            q = _nearest_float(mu[k, j])
            if q != 0
                B[k, :] .-= q .* B[j, :]
                U[k, :] .-= q .* U[j, :]
                for l in 1:j
                    mu[k, l] -= q * mu[j, l]
                end
            end
        end
        if norms[k] >= (d - mu[k, k - 1]^2) * norms[k - 1]
            k += 1
        else
            B[[k - 1, k], :] = B[[k, k - 1], :]
            U[[k - 1, k], :] = U[[k, k - 1], :]
            mu, norms = _gso_float(B)
            k = max(k - 1, 2)
        end
    end
    B, U
end

function _lll_verify(A, B, U, d)
    try
        U * A == B || return false
        abs(det_exact(U)) == 1 || return false
        mu, norms = gram_schmidt_exact(B)
        for i in 1:size(B, 1), j in 1:(i - 1)
            abs(mu[i, j]) <= BigInt(1) // BigInt(2) || return false
        end
        for k in 2:size(B, 1)
            norms[k] >= (d - mu[k, k - 1]^2) * norms[k - 1] || return false
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
