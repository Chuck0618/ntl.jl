using NTL
import Base: zero, one, +, -, *
primes = [2,3,5,7,11, NTL.zz_p_Max]
NTL.isprime.(primes)


@inline function AddMod0(a::Int,b::Int, n::Int)::Int
    return Int(mod(Int128(a)+Int128(b),n))
 end

function test_zz_p(p::Int)
    total = true; 
    for i = 0:p-1, j = 0:p-1
        # @show (i,j)
       total = (AddMod0(i,j,p) == NTL.AddMod(i,j,p));
       if (total == false)
            return false;
       end
    end
    return total;
end

const Int4size = 2^4;
const Int4Max = 2^3-1;
const Int4Min = -2^3;
mutable struct Int4
    rep::Int
    function Int4(a::Int)
        b= mod(a, Int4size)
        c= b>Int4Max ? b-Int4size : b;
        new(c)
    end
end
function Base.+(a::Int4, b::Int4)
    c= a.rep+b.rep;
    c>Int4Max ? c-Int4size : c
end

## [1:3] and 1:3 are different