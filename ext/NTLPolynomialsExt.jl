module NTLPolynomialsExt
using NTL
import Polynomials
Polynomials.Polynomial(f::ModPoly)=Polynomials.Polynomial(coefficients(f))
function NTL.ModPoly(p::Polynomials.Polynomial;ring=nothing)
    c=Polynomials.coeffs(p)
    T=eltype(c)
    F=if !isnothing(ring)
        ring
    elseif T<:Integer
        BigInt
    elseif T<:Rational
        Rational{BigInt}
    elseif T<:Union{GF2,PrimeField}
        T
    elseif isempty(c)
        throw(ArgumentError("supply ring for a zero polynomial with runtime parent"))
    else
        parent(first(c))
    end
    poly(F,c)
end

# NTL operations on external polynomials preserve the external variable.
external(f,p)=Polynomials.Polynomial(coefficients(f),Polynomials.indeterminate(p))
function samevar(p,q)
    Polynomials.indeterminate(p)==Polynomials.indeterminate(q) || throw(ArgumentError("polynomial variables differ"))
end
function NTL.DivRem(p::Polynomials.Polynomial,q::Polynomials.Polynomial)
    samevar(p,q)
    a,b=divrem(ModPoly(p),ModPoly(q))
    external(a,p),external(b,p)
end
function NTL.GCD(p::Polynomials.Polynomial,q::Polynomials.Polynomial)
    samevar(p,q)
    external(gcd(ModPoly(p),ModPoly(q)),p)
end
function NTL.resultant(p::Polynomials.Polynomial,q::Polynomials.Polynomial)
    samevar(p,q)
    resultant(ModPoly(p),ModPoly(q))
end
function NTL.subresultants(p::Polynomials.Polynomial,q::Polynomials.Polynomial)
    samevar(p,q)
    [external(g,p) for g in subresultants(ModPoly(p),ModPoly(q))]
end
NTL.factor(p::Polynomials.Polynomial)=factor(ModPoly(p))
NTL.evaluate(p::Polynomials.Polynomial,a)=evaluate(ModPoly(p),a)
NTL.derivative(p::Polynomials.Polynomial)=external(derivative(ModPoly(p)),p)
NTL.integrate(p::Polynomials.Polynomial)=external(integrate(ModPoly(p)),p)
NTL.roots(p::Polynomials.Polynomial)=roots(ModPoly(p))
NTL.discriminant(p::Polynomials.Polynomial)=discriminant(ModPoly(p))
Polynomials.derivative(p::ModPoly)=NTL.derivative(p)
Polynomials.integrate(p::ModPoly)=NTL.integrate(p)
Polynomials.roots(p::ModPoly)=NTL.roots(p)
end
