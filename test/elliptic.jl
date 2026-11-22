@testset "Elliptic-curve point counting" begin
    elliptic(x, y) = y^2 - x^3 - 2x - 3

    function elements(F::Type{PrimeField{P}}) where {P}
        F
    end
    function elements(E::ExtensionField)
        F = base_field(E)
        d = extension_degree(E)
        (E(collect(c)) for c in Iterators.product(ntuple(_ -> F, d)...))
    end

    function count_fast(F)
        q = cardinality(F)
        affine = 0
        for x in elements(F)
            rhs = x^3 + 2x + 3
            if iszero(rhs)
                affine += 1
            elseif rhs^div(q - 1, 2) == one(rhs)
                affine += 2
            end
        end
        affine + 1
    end

    function count_brute(F)
        affine = 0
        for x in elements(F), y in elements(F)
            iszero(elliptic(x, y)) && (affine += 1)
        end
        affine + 1
    end

    for k in 1:4
        F = k == 1 ? PrimeField(5) : extension_field(PrimeField(5), BuildIrred(PrimeField(5), k; rng = MersenneTwister(k)))
        fast = count_fast(F)
        brute = count_brute(F)
        @test fast == brute
        expected = (7, 25, 127, 625)[k]
        @test fast == expected
    end
end
