using NTL

F = PrimeField(97)             # Validate the prime and obtain the element type.
a = F(120)               # Reduced to 23 modulo 97.
b = F(-8)                # Reduced to 89 modulo 97.

println("a = ", a)
println("b = ", b)
println("a + b = ", a + b)
println("a - b = ", a - b)
println("a * b = ", a * b)
println("a / b = ", a / b)
println("a^100 = ", a^100)
println("a^(-1) = ", a^-1)

@assert a * inv(a) == one(F)
@assert a^0 == one(F)
@assert rep(a) == 23
@assert modulus(a) == 97

# Mixed integer arithmetic reduces the integer automatically.
@assert a + 100 == F(220)
@assert 2a == a + a

# Low-level routines use reduced Int representatives.
@assert AddMod(96, 2, 97) == 1
@assert MulMod(20, 30, 97) == 18
@assert InvMod(20, 97) == 34
@assert PowerMod(20, -1, 97) == 34

# A machine-word prime field is iterable in canonical representative order.
G = PrimeField(5)
@assert length(G) == 5
@assert collect(G) == G.(0:4)
for (i, x) in enumerate(G)
    @assert rep(x) == i - 1
end
println("PrimeField{5} iterates as ", collect(G))

println("All basic arithmetic checks passed.")
