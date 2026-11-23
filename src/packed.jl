export ZZVec, GF2XVec, packedfieldvec, packed_dot
"""Fixed-length integer/polynomial vector with consecutive, little-endian UInt64 limbs per slot."""
struct LimbVector{T} <: AbstractVector{T}
    words::Vector{UInt64}
    stride::Int
    function LimbVector{T}(v::AbstractVector; limbs=nothing) where {T}
        T in (BigInt,GF2X) || throw(ArgumentError("unsupported coefficient type"))
        raw(x)=T===BigInt ? big(x) : bits(GF2X(x))
        width=maximum((ndigits(abs(raw(x)),base=2)+(T===BigInt) for x in v);init=1)
        stride=isnothing(limbs) ? cld(width,64) : Int(limbs)
        stride>0 || throw(ArgumentError("positive limb count required"))
        a=new{T}(zeros(UInt64,Base.checked_mul(length(v),stride)),stride)
        for (i,x) in enumerate(v)
            a[i]=x
        end
        a
    end
end
const ZZVec=LimbVector{BigInt}
const GF2XVec=LimbVector{GF2X}
GF2X(x::GF2X)=x
Base.size(a::LimbVector)=(length(a.words)÷a.stride,)
Base.IndexStyle(::Type{<:LimbVector})=IndexLinear()
function Base.getindex(a::LimbVector{T},i::Int) where T
    checkbounds(a,i)
    last=i*a.stride
    x=big(0)
    for j in last:-1:last-a.stride+1
        x=(x<<64)+a.words[j]
    end
    if T===BigInt && a.words[last]>>63==1
        x-=big(1)<<(64*a.stride)
    end
    T(x)
end
function Base.setindex!(a::LimbVector{T},value,i::Int) where T
    checkbounds(a,i)
    x=T===BigInt ? BigInt(value) : bits(GF2X(value))
    cap=big(1)<<(64*a.stride-(T===BigInt))
    (T===BigInt ? -cap<=x<cap : 0<=x<cap) || throw(InexactError(:setindex!,T,x))
    for j in (i-1)*a.stride+1:i*a.stride
        a.words[j]=x%UInt64
        x >>= 64
    end
    a
end
Base.copy(a::LimbVector{T}) where T=LimbVector{T}(collect(a);limbs=a.stride)
Base.getindex(a::LimbVector{T},r::AbstractUnitRange{Int}) where T=LimbVector{T}([a[i] for i in r];limbs=a.stride)
function packed_dot(a::LimbVector{T},b::LimbVector{T}) where T
    length(a)==length(b) || throw(DimensionMismatch())
    sum((a[i]*b[i] for i in eachindex(a));init=T(0))
end

struct PackedResidues{T,E,F} <: AbstractVector{E}
    ring::F
    data::Vector{T}
end
Base.size(v::PackedResidues)=size(v.data)
Base.IndexStyle(::Type{<:PackedResidues})=IndexLinear()
Base.getindex(v::PackedResidues,i::Int)=v.ring(v.data[i])
@inline Base.getindex(v::PackedResidues{T,PrimeField{P}},i::Int) where {T,P} = _from_integer(Val(P),v.data[i])
function Base.setindex!(v::PackedResidues,x,i::Int)
    v.data[i]=rep(v.ring(x))
    v
end
"""Pack residue arrays, choosing the smallest signed storage that holds p-1."""
function packedfieldvec(F,values;storage=nothing)
    p=modulus(F)
    T=if !isnothing(storage)
        storage
    else
        types=(Int8,Int16,Int32,Int64,Int128)
        i=findfirst(t->p-1<=typemax(t),types)
        isnothing(i) ? BigInt : types[i]
    end
    T<:Integer && (T===BigInt || p-1<=typemax(T)) || throw(ArgumentError("storage cannot represent modulus"))
    data=T[rep(F(x)) for x in values]
    PackedResidues{T,typeof(zero(F)),typeof(F)}(F,data)
end
function +(a::LimbVector{T},b::LimbVector{T}) where T
    length(a)==length(b) || throw(DimensionMismatch())
    LimbVector{T}([a[i]+b[i] for i in eachindex(a)])
end
function -(a::LimbVector{T},b::LimbVector{T}) where T
    length(a)==length(b) || throw(DimensionMismatch())
    LimbVector{T}([a[i]-b[i] for i in eachindex(a)])
end
packedfieldvec(p::Integer,v;kw...) = packedfieldvec(p<=PRIME_FIELD_MAX && ProbPrime(p) ? PrimeField(p) : BigModRing(p),v;kw...)
