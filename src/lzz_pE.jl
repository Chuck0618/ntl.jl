import Base: zero, one, +, -, *, /, //, ^,==, inv, iszero, isone,convert, show,rand
import Primes: isprime
# include("lzz_p.jl")
include("lzz_pX.jl")

struct zz_p_prime{T}
    _rep::zz_pX{T}
    
end