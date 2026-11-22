using NTL

"""Solve a square nonsingular system over F_p by Gaussian elimination."""
function solve_field(A::Matrix{PrimeField{P}}, b::Vector{PrimeField{P}}) where {P}
    n = length(b)
    size(A) == (n, n) || throw(DimensionMismatch("A must be square and match b"))
    aug = hcat(A, b)

    for col in 1:n
        pivot = findfirst(row -> !iszero(aug[row, col]), col:n)
        pivot === nothing && throw(ArgumentError("the matrix is singular"))
        pivot += col - 1
        if pivot != col
            aug[[col, pivot], :] = aug[[pivot, col], :]
        end

        scale = inv(aug[col, col])
        aug[col, :] .*= scale
        for row in 1:n
            row == col && continue
            factor = aug[row, col]
            if !iszero(factor)
                aug[row, :] .-= factor .* aug[col, :]
            end
        end
    end
    aug[:, n + 1]
end

F = PrimeField(101)
A = F.([2 3 1; 4 1 5; 3 7 2])
b = F.([1, 2, 3])
x = solve_field(A, b)

println("solution x = ", x)
println("A*x = ", A * x)
println("b = ", b)
@assert A * x == b
println("The linear-system example passed.")
