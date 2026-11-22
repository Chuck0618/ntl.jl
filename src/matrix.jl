# Generic exact matrix algorithms ----------------------------------------------

export det, matrix_inverse, solve, rank, nullspace, rref, det_exact

function _square_check(A::AbstractMatrix)
    m, n = size(A)
    m == n || throw(DimensionMismatch("the matrix must be square"))
    m > 0 || throw(ArgumentError("an empty matrix has no prototype field element"))
    n
end

"""Determinant over any supported field, using pivoting and Gaussian elimination."""
function _det_generic(A::AbstractMatrix, n::Int)
    M = Matrix(A)
    sign = 1
    factors = [zero(A[1, 1]) for _ in 1:n]
    for k in 1:n
        pivot = findfirst(i -> !iszero(M[i, k]), k:n)
        pivot === nothing && return zero(M[1, 1])
        pivot += k - 1
        if pivot != k
            @inbounds for j in 1:n
                M[k, j], M[pivot, j] = M[pivot, j], M[k, j]
            end
            sign = -sign
        end
        # The same pivot inverse is shared by every row below the pivot.
        pivot_inv = inv(M[k, k])
        @inbounds for i in (k + 1):n
            factors[i] = M[i, k] * pivot_inv
        end
        # Julia matrices are column-major: keep the column fixed in the inner loop.
        @inbounds for j in (k + 1):n
            pivot_row = M[k, j]
            for i in (k + 1):n
                factor = factors[i]
                !iszero(factor) && (M[i, j] -= factor * pivot_row)
            end
        end
    end
    result = M[1, 1]
    for k in 2:n
        result *= M[k, k]
    end
    sign == 1 ? result : -result
end

# After a refresh, w elimination steps require w*(P-1)^2+(P-1)
# to fit in Int. Divide twice to compute the bound without overflow.
_det_refresh_interval(P::Int) = (typemax(Int) - (P - 1)) ÷ (P - 1) ÷ (P - 1)

function _det_primefield_delayed(A::AbstractMatrix{PrimeField{P}}, n::Int, interval::Int) where {P}
    M = Matrix{Int}(undef, n, n)
    @inbounds for (i, a) in enumerate(A)
        M[i] = rep(a)
    end
    sign = 1
    factors = Vector{Int}(undef, n)
    for k in 1:n
        pivot = 0
        @inbounds for i in k:n
            if !iszero(mod(M[i, k], P))
                pivot = i
                break
            end
        end
        pivot == 0 && return _raw(Val{P}(), 0)
        if pivot != k
            @inbounds for j in 1:n
                M[k, j], M[pivot, j] = M[pivot, j], M[k, j]
            end
            sign = -sign
        end
        pivot_inv = _invmod_reduced(mod(M[k, k], P), P)
        @inbounds for i in (k + 1):n
            factors[i] = _mulmod_field(mod(M[i, k], P), pivot_inv, Val{P}())
        end
        if k < n && k % interval == 0
            @inbounds for j in (k + 1):n
                for i in (k + 1):n
                    M[i, j] = mod(M[i, j], P)
                end
            end
        end
        @inbounds for j in (k + 1):n
            pivot_row = mod(M[k, j], P)
            # Dense matrices rarely benefit from skipping zero factors; keeping
            # contiguous row indices lets LLVM vectorize the integer update.
            @simd for i in (k + 1):n
                M[i, j] -= factors[i] * pivot_row
            end
        end
    end
    result = 1
    @inbounds for k in 1:n
        result = _mulmod_field(result, mod(M[k, k], P), Val{P}())
    end
    _raw(Val{P}(), sign == 1 ? result : (iszero(result) ? 0 : P - result))
end

function det(A::AbstractMatrix{PrimeField{P}}) where {P}
    n = _square_check(A)
    interval = _det_refresh_interval(P)
    interval >= 2 && return _det_primefield_delayed(A, n, min(n, interval))
    _det_generic(A, n)
end

function det(A::AbstractMatrix)
    _det_generic(A, _square_check(A))
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
    # Give each cell unique storage before using in-place GMP operations.
    M = [Base.GMP.MPZ.set!(BigInt(0), BigInt(A[i, j])) for i in axes(A, 1), j in axes(A, 2)]
    n == 1 && return M[1, 1]
    sign = 1
    previous = BigInt(1)
    left, right = BigInt(0), BigInt(0)
    for k in 1:(n - 1)
        if iszero(M[k, k])
            pivot = findfirst(i -> !iszero(M[i, k]), (k + 1):n)
            pivot === nothing && return BigInt(0)
            pivot += k
            for j in k:n
                M[k, j], M[pivot, j] = M[pivot, j], M[k, j]
            end
            sign = -sign
        end
        pivot_value = M[k, k]
        for j in (k + 1):n, i in (k + 1):n
            _bareiss_update!(M[i, j], pivot_value, M[i, k], M[k, j], previous, left, right)
        end
        previous = pivot_value
    end
    sign * M[n, n]
end

@inline function _bareiss_update!(cell, pivot, column, row, previous, left, right)
    Base.GMP.MPZ.mul!(left, pivot, cell)
    Base.GMP.MPZ.mul!(right, column, row)
    Base.GMP.MPZ.sub!(left, left, right)
    Base.GMP.MPZ.tdiv_q!(cell, left, previous)
end
