import Base: zero, one, +, -, *, /, //, ^,==, inv, iszero, isone,convert, show,rand
import Primes: isprime

struct NTL_INIT_zz_p end # type for init
const _init_zz_p_nocheck=NTL_INIT_zz_p()
const zz_p_Max=begin# maximal prime number for package
    i = typemax(Int) 
    while(true)
    global i;
    if(isprime(i))
        break; 
    end;
    i = i-1;
    end
    i
end
const zz_p_Min=2; # minimal prime number, clearly 2

# 备用素数type 
struct zz_p_prime 
    _prime::Int
    function zz_p_prime(a::Int)
        new(a)
    end
end
function convert(a::Int,::zz_p_prime)
    # test prime ??
    return zz_p_prime(a)
end


"""
struct zz_p{T} is the type of finite fields

    
    ## Examples:
    julia> zz_p{7} # The type of finite fields of cardinality 7
    julia> zz_p{3}(4) # Element 4 in zz_p{3} which equals 1 in zz_p{3}
    julia> zz_p{3}() # Element 0 in zz_p{3}
"""
mutable struct zz_p{T}
    _rep::Int 
    function zz_p{T}() where {T} 
        new(0::Int)
    end
    # If we know a is representable
    function zz_p{T}(a::Int, ::NTL_INIT_zz_p) where {T} 
        new(a)
    end
    function zz_p{T}(a::Int) where {T} 
        new(mod(a, T))
    end
end



## template
"""
Generate the type zz_p{T} of prime fields

    ## Examples:
    julia> zz_p(3)
    julia> zz_p{3} # they are equivalent
"""
function zz_p(T::Int)
    if T > 1 && isprime(T)
    return zz_p{T}
    end
    error("The parameter is NOT a prime!")
end

# ***************************************************************
#                          Information and Conversion
# ***************************************************************
## info of type zz_p
function info(::Type{zz_p})
    return Dict(
        "Discript"=>"This module is a prime field with small p", "zz_p_Max"=>zz_p_Max,
        "zz_p_Min"=>zz_p_Min)
end
## info of type zz_p{T}
function modulus(::Type{zz_p{T}}) where {T}
    return T
end
function zero(::Type{zz_p{T}}) where {T}
    return zz_p{T}(0)
end
function one(::Type{zz_p{T}}) where {T}
    return zz_p{T}(1, _init_zz_p_nocheck)
end
function info(::Type{zz_p{T}}) where {T}
    return Dict("modulus"=>modulus(zz_p{T}), "zero"=>zero(zz_p{T}), "one"=>one(zz_p{T}))
end
## info of instance of zz_p{T}
function info(a::zz_p{T}) where {T}
    return Dict("reprentation"=>a._rep, "field"=>"zz_"*string(T))
end

## convert
function convert(::Type{zz_p{T}}, a::Int) where {T}
    return zz_p{T}(a)
end
function convert(::Type{Int}, a::zz_p{T}) where {T}
    return a._rep
end
## swap
@inline function swap!(x::zz_p{T}, y::zz_p{T}) where {T}
    x,y = y,x ;
    nothing
end

## inline funcitons
#const NTL_BITS_PER_LONG=Sys.WORD_SIZE
# @inline function cast_unsigned(a::Int)
#     return Core.bitcast(UInt, a)
# end
# @inline function sp_SignMask(a::Int)
#    return -(cast_unsigned(a) >> (NTL_BITS_PER_LONG-1));
# end
@inline function sp_CorrectExcess(a::Int, n::Int)
   return a-n >= 0 ? a-n : a
end

## add
# @inline function AddMod0(a::Int,b::Int, n::Int)::Int
#     return Int(mod(Int128(a)+Int128(b),n))
#  end
@inline function AddMod(a::Int,b::Int, n::Int)::Int
    r = a+b;
    return sp_CorrectExcess(r, n)
 end
# function add!(z::zz_p{T} , x::zz_p{T}, y::zz_p{T}) where {T}
#     z._rep = AddMod(x._rep,y._rep, T)
#     return z
# end
function add(x::zz_p{T}, y::zz_p{T}) where {T}
    Z=AddMod(x._rep,y._rep, T)
    return zz_p{T}(Z, _init_zz_p_nocheck) # since we don't need to check whether Z is within 0 .. p-1
end
+(x::zz_p{T}, y::zz_p{T}) where {T} = add(x,y) 
+(x::zz_p{T}, Y::Int) where {T} = add(x,convert(zz_p{T},Y))
+(X::Int, y::zz_p{T}) where {T} = add(convert(zz_p{T},X),y)

# ***************************************************************
#                          Subtraction
# ***************************************************************
@inline function sp_CorrectDeficit(a::Int, n::Int)
    return a >= 0 ? a : a+n
 end

@inline function SubMod0(a::Int,b::Int, n::Int)::Int
    return mod(a-b, n)
 end
@inline function SubMod(a::Int,b::Int, n::Int)::Int
    r = a-b;
    return sp_CorrectDeficit(r, n)
 end
function sub!(z::zz_p{T} , x::zz_p{T}, y::zz_p{T}) where {T}
    z._rep = SubMod(x._rep,y._rep, T)
    return z
end
function sub(x::zz_p{T}, y::zz_p{T}) where {T}
    Z=SubMod(x._rep,y._rep, T)
    return zz_p{T}(Z, _init_zz_p_nocheck)
end
-(x::zz_p{T}, y::zz_p{T}) where {T} = sub(x,y)
-(x::zz_p{T}, Y::Int) where {T} = sub(x,convert(zz_p{T},Y))
-(X::Int, y::zz_p{T}) where {T} = sub(convert(zz_p{T},X),y)

## negate
# x = -a
@inline function NegateMod0(a::Int, n::Int)::Int
    return mod((0-a) , n)
end

@inline function NegateMod(a::Int, n::Int)::Int
    # mod((0- a) , n)
      return n-a;
end
function negate!(x::zz_p{T}, a::zz_p{T}) where {T}
    x._rep=NegateMod(a._rep,T)
    return x
end
function negate(a::zz_p{T}) where {T} 
    Z = NegateMod(a._rep,T)
    return zz_p{T}(Z, _init_zz_p_nocheck)
end
-(a::zz_p{T}) where {T} = negate(a)

# ***************************************************************
#                          Multiplication
# ***************************************************************
@inline function MulMod(a::Int, b::Int, n::Int)::Int
    Int(mod(Int128(a)*Int128(b),n))
end
function mul!(z::zz_p{T} , x::zz_p{T}, y::zz_p{T}) where {T}
    z._rep = MulMod(x._rep,y._rep, T)
    return z
end
function mul(x::zz_p{T}, y::zz_p{T}) where {T}
    Z=MulMod(x._rep,y._rep, T)
    return zz_p{T}(Z, _init_zz_p_nocheck)
end

*(x::zz_p{T}, y::zz_p{T}) where {T} = mul(x,y)
*(x::zz_p{T}, Y::Int) where {T} = mul(x, convert(zz_p{T},Y))
*(X::Int, y::zz_p{T}) where {T} = mul(convert(zz_p{T},X), y)


# ***************************************************************
#                          Inverse
# ***************************************************************
## inv 
# 1/a mod n
function InvMod(a::Int, n::Int)::Int
    d,u,v = gcdx(a, n)
    if d!= 1
        error("InvMod: inverse undefined")
    elseif u< 0
        return u+n
    end
    return u
end
function inv!(z::zz_p{T}, x::zz_p{T}) where {T}
    z._rep=InvMod(x._rep, T)
    return z
end
function inv(x::zz_p{T}) where {T}
    Z = InvMod(x._rep, T)
    return zz_p{T}(Z, _init_zz_p_nocheck)
end

##div
function div!(z::zz_p{T} , x::zz_p{T}, y::zz_p{T}) where {T}
    z._rep = MulMod(x._rep, InvMod(y._rep, T), T)
    return z
end
function div(x::zz_p{T}, y::zz_p{T}) where {T}
    Z=MulMod(x._rep, InvMod(y._rep, T), T)
    return zz_p{T}(Z, _init_zz_p_nocheck)
end
/(x::zz_p{T}, y::zz_p{T}) where {T} = div(x,y)
/(x::zz_p{T}, Y::Int) where {T} = div(x, convert(zz_p{T},Y))
/(X::Int, y::zz_p{T}) where {T} = div(convert(zz_p{T},X), y)


#***************************************************************
#                          comparison
#***************************************************************

@inline function iszero(x::zz_p{T}) where {T}
    return x._rep == 0
end
@inline function isone(x::zz_p{T}) where {T}
    return x._rep == 1
end
@inline function ==(x::zz_p{T}, y::zz_p{T}) where {T}
    return x._rep == y._rep
end
@inline function ==(x::zz_p{T}, Y::Int) where {T}
    return x._rep == zz_p{T}(Y)
end
@inline function ==(X::Int, y::zz_p{T}) where {T}
    return zz_p{T}(X) == y._rep
end
## random numbers
function rand(::Type{zz_p{T}}) where {T}
    return zz_p{T}(rand(Int))
end
function rand(::Type{zz_p{T}}, k::Int) where {T}
    return zz_p{T}.(rand(Int,k))
end


#***************************************************************
#                          power
#***************************************************************
## power x = a^e in zz_p
function PowerMod(a::Int, e::Int , n::Int)
    if e < 0
        return InvMod(PowerMod(a,-e,n),n)
    end
        x::Int = 1;
        y::Int = a;
        while (e != 0)
            if (e & 1) != 0 
                    x = MulMod(x, y, n)
            end
            y = MulMod(y, y, n);
            e = e >> 1;
        end
    return x;
end

@inline function power!(z::zz_p{T}, x::zz_p{T}, e::Int) where {T}
    z._rep= PowerMod(x._rep, e, T)
    return z
end

@inline function power(x::zz_p{T}, e::Int) where {T}
    Z= PowerMod(x._rep, e, T)
    return zz_p{T}(Z, _init_zz_p_nocheck)
end

^(x::zz_p{T}, e::Int) where {T} = power(x,e)

#***************************************************************
#                          iostream
#***************************************************************
## show
show(io::IO, x::zz_p{T}) where {T} = print(io, "zz_", T, "(", x._rep ,")" )
show(io::IO, ::Type{zz_p{T}}) where {T} = print(io, "Finite Field: zz_", modulus(zz_p{T}))


