export RR, xdouble, quad_float, mat_RR, isolate_real_roots, BasicThreadPool, parallel_map, Lazy, LazyTable
const RR=BigFloat
const mat_RR=Matrix{BigFloat}
"""Extended-exponent software real rounded to P significant binary digits."""
struct SoftReal{P} <: Real
    value::BigFloat
    SoftReal{P}(x::Real) where P=new{P}(BigFloat(x;precision=P))
end
const xdouble=SoftReal{53}
const quad_float=SoftReal{106}
Base.BigFloat(x::SoftReal;precision::Integer=Base.precision(BigFloat))=BigFloat(x.value;precision)
Base.Float64(x::SoftReal)=Float64(x.value)
Base.convert(::Type{T},x::T) where {T<:SoftReal}=x
Base.convert(::Type{T},x::Real) where {T<:SoftReal}=T(x)
Base.promote_rule(::Type{T},::Type{S}) where {T<:SoftReal,S<:Integer}=T
Base.promote_rule(::Type{SoftReal{P}},::Type{SoftReal{Q}}) where {P,Q}=SoftReal{max(P,Q)}
Base.show(io::IO,x::SoftReal)=show(io,x.value)
Base.:(==)(x::SoftReal,y::SoftReal)=x.value==y.value
Base.isless(x::SoftReal,y::SoftReal)=isless(x.value,y.value)
Base.:(<)(x::SoftReal,y::SoftReal)=x.value<y.value
Base.:(<=)(x::SoftReal,y::SoftReal)=x.value<=y.value
Base.hash(x::SoftReal,h::UInt)=hash(x.value,h)
# MPFR writes into a freshly allocated destination of the requested precision.
# The global BigFloat precision/rounding context is never changed here.
for (op,name) in ((:+,:mpfr_add),(:-,:mpfr_sub),(:*,:mpfr_mul),(:/,:mpfr_div))
    @eval function Base.$op(x::SoftReal{P},y::SoftReal{P}) where P
        z=BigFloat(0;precision=P)
        ccall(($(QuoteNode(name)),Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Cint),z,x.value,y.value,0)
        SoftReal{P}(z)
    end
end
for (op,name) in ((:-,:mpfr_neg),(:sqrt,:mpfr_sqrt),(:exp,:mpfr_exp),(:log,:mpfr_log),(:sin,:mpfr_sin),(:cos,:mpfr_cos),(:tan,:mpfr_tan))
    @eval function Base.$op(x::SoftReal{P}) where P
        z=BigFloat(0;precision=P)
        ccall(($(QuoteNode(name)),Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Cint),z,x.value,0)
        SoftReal{P}(z)
    end
end

"""Disjoint open rational intervals, one per distinct real root (Sturm)."""
function isolate_real_roots(f::ModPoly{BigInt};width=1//1024)
    width>0 || throw(DomainError(width))
    iszero(f) && throw(DomainError(f))
    Q=Rational{BigInt}
    out=Tuple{Q,Q}[]
    degree(f)==0 && return out
    p=poly(Q,coefficients(f))
    p=div(p,gcd(p,derivative(p)))
    seq=[p,derivative(p)]
    while true
        r=-rem(seq[end-1],seq[end])
        iszero(r) && break
        push!(seq,r)
    end
    function changes(x)
        s=filter(!iszero,[sign(evaluate(q,x)) for q in seq])
        count(i->s[i]!=s[i-1],2:length(s))
    end
    limit=Q(2+ceil(BigInt,maximum(abs,p._coeffs)/abs(LeadCoeff(p))))
    stack=[(-limit,limit)]
    while !isempty(stack)
        a,b=pop!(stack)
        n=changes(a)-changes(b)
        n==0 && continue
        if n==1 && b-a<=width
            push!(out,(a,b))
            continue
        end
        mid=(a+b)/2
        while iszero(evaluate(p,mid))
            mid=(a+mid)/2
        end
        push!(stack,(mid,b),(a,mid))
    end
    out
end

struct BasicThreadPool
    workers::Int
    function BasicThreadPool(workers::Integer=Threads.nthreads())
        workers>0 || throw(DomainError(workers))
        new(Int(workers))
    end
end
"""Bound concurrent jobs while preserving input order; propagate task errors."""
function parallel_map(f,pool::BasicThreadPool,items)
    xs=collect(items)
    isempty(xs) && return map(f,xs)
    workers=min(pool.workers,length(xs))
    tasks=[Threads.@spawn [f(xs[i]) for i in k:workers:length(xs)] for k in 1:workers]
    # Join every child before reporting an error.
    @sync for t in tasks
        @async wait(t)
    end
    batches=fetch.(tasks)
    [batches[mod1(i,workers)][cld(i,workers)] for i in eachindex(xs)]
end
struct LazyTable{K,V}
    values::Dict{K,V}
    lock::ReentrantLock
end
LazyTable{K,V}() where {K,V}=LazyTable(Dict{K,V}(),ReentrantLock())
Base.get!(f::Function,t::LazyTable,key)=lock(()->get!(f,t.values,key),t.lock)
struct Lazy{F,V}
    builder::F
    table::LazyTable{Nothing,V}
end
Lazy(f::F,::Type{V}=Any) where {F,V}=Lazy(f,LazyTable{Nothing,V}())
Base.getindex(l::Lazy)=get!(l.builder,l.table,nothing)

# User-owned plans bound cache lifetime to the batch using them.
export ArithmeticCache, cached_field, cached_ntt_convolve, cached_zech
struct ArithmeticCache
    fields::LazyTable{Int,DataType}
    roots::LazyTable{Tuple{Int,Int},Int}
    twiddles::LazyTable{Tuple{Int,Int,Int},Vector{Int}}
    zech::LazyTable{Any,Any}
end
ArithmeticCache()=ArithmeticCache(LazyTable{Int,DataType}(),LazyTable{Tuple{Int,Int},Int}(),
    LazyTable{Tuple{Int,Int,Int},Vector{Int}}(),LazyTable{Any,Any}())
cached_field(c::ArithmeticCache,p::Int)=get!(() -> PrimeField(p),c.fields,p)
function cached_zech(c::ArithmeticCache,E::ExtensionField)
    key=(typeof(base_field(E)),characteristic(E),Tuple(rep.(coefficients(modulus(E)))),E.symbol)
    get!(() -> enable_zech_multiplication(E),c.zech,key)
end
function _plan_root(c::ArithmeticCache,v::Val{P},n) where P
    _raw(v,get!(() -> rep(_convolution_root(v,n)), c.roots,(P,n)))
end
function _plan_twiddles(c::ArithmeticCache,w::PrimeField{P},n) where P
    get!(() -> _convolution_twiddles(w,n), c.twiddles,(P,n,rep(w)))
end
cached_ntt_convolve(c::ArithmeticCache,a::Vector{PrimeField{P}},b::Vector{PrimeField{P}}) where P = _ntt_convolve_plan(c,a,b)
