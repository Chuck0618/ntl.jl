# Exact LLL lattice reduction --------------------------------------------------

export gram_schmidt_exact, lll

function _gram_bareiss(B::AbstractMatrix{<:Integer})
    n, m = size(B)
    n > 0 && m > 0 || throw(ArgumentError("the basis must be nonempty"))
    Z = Matrix{BigInt}(B)
    gram = Matrix{BigInt}(undef, n, n)
    left, right = BigInt(0), BigInt(0)
    for j in 1:n, i in j:n
        sum = BigInt(0)
        for k in 1:m
            Base.GMP.MPZ.mul!(left, Z[i, k], Z[j, k])
            Base.GMP.MPZ.add!(sum, sum, left)
        end
        gram[i, j] = sum
    end
    previous = BigInt(1)
    # Bareiss elimination on the lower triangle of the Gram matrix.
    # The pivot is the leading principal determinant, and each entry
    # below it is that pivot times the corresponding GSO coefficient.
    for j in 1:n
        pivot = gram[j, j]
        iszero(pivot) && throw(ArgumentError("the input rows must be linearly independent"))
        for k in (j + 1):n, i in k:n
            _bareiss_update!(gram[i, k], pivot, gram[i, j], gram[k, j], previous, left, right)
        end
        previous = pivot
    end
    gram
end

"""Exact rational Gram–Schmidt data for a full-rank row basis."""
function gram_schmidt_exact(B::AbstractMatrix{<:Integer})
    gram = _gram_bareiss(B)
    n = size(gram, 1)
    Q = Rational{BigInt}
    mu = zeros(Q, n, n)
    norms = zeros(Q, n)
    previous = BigInt(1)
    for j in 1:n
        pivot = gram[j, j]
        norms[j] = pivot // previous
        for i in (j + 1):n
            mu[i, j] = gram[i, j] // pivot
        end
        previous = pivot
    end
    mu, norms
end

function _nearest_integer(x::Rational{BigInt})
    floor(BigInt, x + BigInt(1) // BigInt(2))
end

function _lll_delta(delta)
    d = delta isa Union{Rational,Integer} ? Rational{BigInt}(delta) : rationalize(BigInt, delta)
    BigInt(1) // BigInt(4) < d < 1 || throw(DomainError(delta, "LLL delta must lie strictly between 1/4 and 1"))
    d
end

"""
    lll(A; delta=3//4) -> (basis, transform)

Exact-arithmetic LLL for a full-rank integer row basis. Returns a named tuple
whose `basis` and `transform` satisfy `transform * A == basis`.
"""
function _lll_exact_legacy(A::AbstractMatrix{<:Integer}; delta = BigInt(3) // BigInt(4), maxiter::Integer = 100000)
    d = _lll_delta(delta)
    n, m = size(A)
    n > 0 && m > 0 || throw(ArgumentError("the basis must be nonempty"))
    maxiter > 0 || throw(DomainError(maxiter, "maxiter must be positive"))
    B = BigInt.(A)
    U = zeros(BigInt, n, n)
    for i in 1:n
        U[i, i] = 1
    end
    n == 1 && return (basis = B, transform = U)
    mu, norms = gram_schmidt_exact(B)
    k = 2
    iterations = 0
    while k <= n
        iterations += 1
        iterations <= maxiter || error("LLL iteration limit exceeded")
        for j in (k - 1):-1:1
            q = _nearest_integer(mu[k, j])
            if !iszero(q)
                B[k, :] .-= q .* B[j, :]
                U[k, :] .-= q .* U[j, :]
                mu, norms = gram_schmidt_exact(B)
            end
        end
        if norms[k] >= (d - mu[k, k - 1]^2) * norms[k - 1]
            k += 1
        else
            B[[k - 1, k], :] = B[[k, k - 1], :]
            U[[k - 1, k], :] = U[[k, k - 1], :]
            mu, norms = gram_schmidt_exact(B)
            k = max(k - 1, 2)
        end
    end
    (basis = B, transform = U)
end
