# Generic exact matrix algorithms ----------------------------------------------

export det, matrix_inverse, solve, rank, nullspace, rref, det_exact

function _square_check(A::AbstractMatrix)
    m, n = size(A)
    m == n || throw(DimensionMismatch("the matrix must be square"))
    m > 0 || throw(ArgumentError("an empty matrix has no prototype field element"))
    n
end

"""Determinant over any supported field, using pivoting and Gaussian elimination."""
function det(A::AbstractMatrix)
    n = _square_check(A)
    M = Matrix(A)
    sign = 1
    for k in 1:n
        pivot = findfirst(i -> !iszero(M[i, k]), k:n)
        pivot === nothing && return zero(M[1, 1])
        pivot += k - 1
        if pivot != k
            M[[k, pivot], :] = M[[pivot, k], :]
            sign = -sign
        end
        pivot_value = M[k, k]
        for i in (k + 1):n
            factor = M[i, k] / pivot_value
            iszero(factor) && continue
            for j in (k + 1):n
                M[i, j] -= factor * M[k, j]
            end
        end
    end
    result = M[1, 1]
    for k in 2:n
        result *= M[k, k]
    end
    sign == 1 ? result : -result
end

"""Reduced row-echelon form and pivot column indices."""
function rref(A::AbstractMatrix)
    m, n = size(A)
    (m > 0 && n > 0) || throw(ArgumentError("rref needs a nonempty matrix"))
    M = Matrix(A)
    pivots = Int[]
    row = 1
    for col in 1:n
        row > m && break
        pivot = findfirst(i -> !iszero(M[i, col]), row:m)
        pivot === nothing && continue
        pivot += row - 1
        pivot != row && (M[[row, pivot], :] = M[[pivot, row], :])
        M[row, :] .*= inv(M[row, col])
        for i in 1:m
            i == row && continue
            factor = M[i, col]
            iszero(factor) && continue
            M[i, :] .-= factor .* M[row, :]
        end
        push!(pivots, col)
        row += 1
    end
    M, pivots
end

"""Rank over an exact field."""
function rank(A::AbstractMatrix)
    _, pivots = rref(A)
    length(pivots)
end

"""Solve `A*x=b` over a field for a vector or multiple matrix right-hand sides."""
function solve(A::AbstractMatrix, b::Union{AbstractVector,AbstractMatrix})
    m, n = size(A)
    m == n || throw(DimensionMismatch("solve expects a square coefficient matrix"))
    size(b, 1) == n || throw(DimensionMismatch("right-hand side has incompatible size"))
    n > 0 || throw(ArgumentError("empty system"))
    aug = hcat(A, b)
    R, pivots = rref(aug)
    pivots == collect(1:n) || throw(ArgumentError("the coefficient matrix is singular"))
    R[:, (n + 1):end]
end

function solve(A::AbstractMatrix, b::AbstractVector)
    m, n = size(A)
    m == n || throw(DimensionMismatch("solve expects a square coefficient matrix"))
    length(b) == n || throw(DimensionMismatch("right-hand side has incompatible length"))
    n > 0 || throw(ArgumentError("empty system"))
    aug = hcat(A, b)
    R, pivots = rref(aug)
    pivots == collect(1:n) || throw(ArgumentError("the coefficient matrix is singular"))
    R[:, n + 1]
end

"""Inverse of a nonsingular square matrix over a field."""
function matrix_inverse(A::AbstractMatrix)
    n = _square_check(A)
    proto = A[1, 1]
    I = [i == j ? one(proto) : zero(proto) for i in 1:n, j in 1:n]
    solve(A, I)
end

"""A vector basis for the right nullspace of `A`."""
function nullspace(A::AbstractMatrix)
    m, n = size(A)
    R, pivots = rref(A)
    free = setdiff(1:n, pivots)
    proto = A[1, 1]
    basis = Vector{Vector{typeof(proto)}}()
    for f in free
        v = [zero(proto) for _ in 1:n]
        v[f] = one(proto)
        for (row, col) in enumerate(pivots)
            v[col] = -R[row, f]
        end
        push!(basis, v)
    end
    basis
end

"""Exact integer determinant by fraction-free Bareiss elimination."""
function det_exact(A::AbstractMatrix{<:Integer})
    n = _square_check(A)
    M = BigInt.(A)
    n == 1 && return M[1, 1]
    sign = 1
    previous = BigInt(1)
    for k in 1:(n - 1)
        if iszero(M[k, k])
            pivot = findfirst(i -> !iszero(M[i, k]), (k + 1):n)
            pivot === nothing && return BigInt(0)
            pivot += k
            M[[k, pivot], :] = M[[pivot, k], :]
            sign = -sign
        end
        pivot_value = M[k, k]
        for i in (k + 1):n, j in (k + 1):n
            numerator = M[i, j] * pivot_value - M[i, k] * M[k, j]
            M[i, j] = div(numerator, previous)
        end
        previous = pivot_value
    end
    sign * M[n, n]
end
