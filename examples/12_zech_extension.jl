using NTL
using Statistics
using Printf

F = PrimeField(3)
E = extension_field(F, poly(F, [2, 1, 1]))
Z = zech_field(E)
x = ZechElem(Z, Z._generator)

@assert cardinality(Z) == 9
@assert x * inv(x) == one(Z)
@assert x^8 == one(Z)
@assert muladd(x, x, one(x)) == x * x + one(x)
@assert isbitstype(typeof(x))

println("Zech element = ", x)
println("inverse = ", inv(x))
println("x^8 = ", x^8)

for fn in (() -> x * x + one(x), () -> muladd(x, x, one(x)))
    fn()
    ts = [@elapsed fn() for _ in 1:8]
    println(1000 * Statistics.median(ts), " ms")
end
