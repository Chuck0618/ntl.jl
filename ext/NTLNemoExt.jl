module NTLNemoExt
using NTL
import Nemo
function factor_zz(f)
    R,_ = Nemo.polynomial_ring(Nemo.ZZ,:x)
    fac = Nemo.factor(R(Nemo.ZZ.(coefficients(f))))
    factors = Tuple{typeof(f),Int}[]
    for (g,e) in fac
        push!(factors,(poly(BigInt,[BigInt(Nemo.coeff(g,i)) for i in 0:Nemo.degree(g)]),e))
    end
    PolyFactorization(BigInt(Nemo.coeff(Nemo.unit(fac),0)),factors)
end
end
