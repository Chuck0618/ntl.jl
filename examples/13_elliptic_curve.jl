using NTL
using Random
using Statistics
using Printf

# Count projective points of y^2 = x^3 + 2x + 3 over F_(5^k).
# The cubic is singular in characteristic five, but the point count is still well-defined.

elliptic(x, y) = y^2 - x^3 - 2x - 3

function make_field(k)
    F0 = PrimeField(5)
    k == 1 && return F0
    extension_field(F0, BuildIrred(F0, k; rng = MersenneTwister(2026 + k)))
end

all_elements(F::Type{PrimeField{P}}) where {P} = F
function all_elements(E::ExtensionField)
    F = base_field(E)
    d = extension_degree(E)
    (E(collect(coeffs)) for coeffs in Iterators.product(ntuple(_ -> F, d)...))
end

function square_target_count(a, q)
    iszero(a) && return 1
    a^div(q - 1, 2) == one(a) ? 2 : 0
end

function count_points_fast(F)
    q = cardinality(F)
    affine = 0
    for x in all_elements(F)
        affine += square_target_count(x^3 + 2x + 3, q)
    end
    (affine = affine, projective = affine + 1)
end

function count_points_bruteforce(F)
    affine = 0
    elements = all_elements(F)
    for x in elements, y in elements
        iszero(elliptic(x, y)) && (affine += 1)
    end
    (affine = affine, projective = affine + 1)
end

function main(max_degree = 6)
    println("E: y^2 = x^3 + 2x + 3 over F_(5^k)")
    println("The projective count includes the point at infinity.")
    for k in 1:max_degree
        F = make_field(k)
        started = time()
        result = count_points_fast(F)
        elapsed = 1000 * (time() - started)
        if cardinality(F) <= 125
            @assert result == count_points_bruteforce(F)
        end
        @printf("k=%d  q=%d  affine=%d  projective=%d  time=%.3f ms\n",
                k, cardinality(F), result.affine, result.projective, elapsed)
    end
    @assert count_points_fast(PrimeField(5)).projective == 7
    println("Elliptic-curve point counting passed.")
end

main()
