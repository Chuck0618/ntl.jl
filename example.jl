using NTL
A=[2 4 6;1 3 3]
H,U=hnf(A)
@assert U*A==H && abs(det_exact(U))==1
s=solve_integer(A,[2,1])
@assert A*s.x==[2,1] && iszero(A*s.kernel)
@show H s.x s.kernel charpoly([1 2;3 4])
E,a=conway_field(2,3)
x=xpoly(E)
f=(x+a)^2*(x+a+1)
z=GF2EXFactoring(f)
@assert reconstruct(z)==f && (x+a,2) in z.factors
@assert minpoly(a)==defining_polynomial(E)
t=xpoly(BigInt)
@assert reconstruct(factor_integer_poly(-6*(2*t+1)^2*(t^2+1)))==-6*(2*t+1)^2*(t^2+1)
A=[28 18 -1 -6;29 21 -7 27;-21 6 11 5;-28 -9 24 -7]
B=lll(A;delta=99//100,mode=:exact).basis
r=bkz(A;block_size=4)
@assert r.verified && r.transform*A==r.basis
println("First-vector squared norm: LLL ",sum(abs2,B[1,:]),", BKZ ",sum(abs2,r.basis[1,:]))
