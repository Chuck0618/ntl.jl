using NTL

F = zzp(101)
A = F.([2 3 1; 4 1 5; 3 7 2])
b = F.([1, 2, 3])
x = solve(A, b)
@assert A * x == b
@assert A * matrix_inverse(A) == F.([1 0 0; 0 1 0; 0 0 1])
println("solution over F_101 = ", x)
println("det(A) = ", det(A))

singular = F.([1 2 3; 2 4 6; 0 1 1])
K = nullspace(singular)
@assert length(K) == 1 && all(iszero, singular * K[1])
println("nullspace basis = ", K)

B0 = BigInt[1 1 1; -1 0 2; 3 5 6]
reduced = lll(B0)
B, U = reduced.basis, reduced.transform
@assert U * B0 == B
@assert abs(det_exact(U)) == 1
println("LLL-reduced basis = ", B)
println("unimodular transform = ", U)
println("Matrix and lattice example passed.")
