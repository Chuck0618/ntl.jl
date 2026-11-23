export charpoly, MinPolyMod, minpoly, hnf, kernel_integer, image_integer, solve_integer

"""Division-free Berkowitz characteristic polynomial, including composite residue rings."""
function charpoly(A::AbstractMatrix; ring=nothing)
    n, m = size(A)
    n == m || throw(DimensionMismatch("expected a square matrix"))
    F = if ring !== nothing
        ring
    elseif eltype(A) <: Integer
        BigInt
    elseif n > 0
        parent(first(A))
    elseif eltype(A) <: Union{PrimeField,GF2} && isconcretetype(eltype(A))
        eltype(A)
    else
        throw(ArgumentError("supply ring for an empty matrix with runtime parent"))
    end
    M = map(F, A)
    c = [one(F)]
    for k in n:-1:1
        # First column of the lower triangular Toeplitz recurrence.
        t = [one(F), -M[k,k]]
        v = M[k+1:n,k]
        for j in k+1:n
            push!(t, -sum((M[k,i]*v[i-k] for i in k+1:n); init=zero(F)))
            j < n && (v = [sum((M[i,l]*v[l-k] for l in k+1:n); init=zero(F)) for i in k+1:n])
        end
        c = [sum((t[i-j+1]*c[j] for j in 1:min(i,length(c))); init=zero(F)) for i in eachindex(t)]
    end
    poly(F, reverse(c))
end

"""First monic dependence among 1,h,h²,… in F[x]/(f), over a field."""
function MinPolyMod(h::ModPoly, f::ModPoly)
    F = _check_poly_parents(h,f)
    degree(f) > 0 || throw(DomainError(f,"positive degree required"))
    F === BigInt && throw(ArgumentError("a coefficient field is required"))
    _validate_factor_parent(F)
    n = degree(f)
    pivots = zeros(Int,n)
    columns = Vector{Vector{typeof(zero(F))}}()
    relations = similar(columns,0)
    v = one(f)
    h = rem(h,f)
    for k in 0:n
        w = [coeff(v,i) for i in 0:n-1]
        relation = fill(zero(F),n+1)
        relation[k+1] = one(F)
        for i in 1:n
            iszero(w[i]) && continue
            j = pivots[i]
            if iszero(j)
                scale = inv(w[i])
                push!(columns,[x*scale for x in w])
                push!(relations,[x*scale for x in relation])
                pivots[i] = length(columns)
                break
            end
            a = w[i]
            w .-= [a*x for x in columns[j]]
            relation .-= [a*x for x in relations[j]]
        end
        all(iszero,w) && return poly(F,relation)
        v = MulMod(v,h,f)
    end
    error("missing quotient dependence")
end
minpoly(a::ExtensionElem) = MinPolyMod(rep(a),modulus(a))
minpoly(a::Union{PrimeField,GF2,ModInt}) = poly(parent(a),[-a,one(a)])

"""Row Hermite form, positive pivots, zero rows last. With transform returns (H,U)."""
function hnf(A::AbstractMatrix{<:Integer}; transform::Bool=true)
    H = BigInt.(A)
    n,m = size(H)
    U = [BigInt(i==j) for i in 1:n,j in 1:(transform ? n : 0)]
    k = 1
    for col in 1:m
        k>n && break
        for i in k+1:n
            b = H[i,col]
            iszero(b) && continue
            a = H[k,col]
            g,s,t = gcdx(a,b)
            for M in (H,U), j in axes(M,2)
                x,y = M[k,j],M[i,j]
                M[k,j],M[i,j] = s*x+t*y, div(a,g)*y-div(b,g)*x
            end
        end
        p = H[k,col]
        iszero(p) && continue
        if p < 0
            H[k,:] .*= -1
            U[k,:] .*= -1
            p = -p
        end
        for i in 1:k-1
            q = fld(H[i,col],p)
            H[i,:] .-= q .* H[k,:]
            U[i,:] .-= q .* U[k,:]
        end
        k += 1
    end
    transform ? (;H,U) : H
end
_nonzero_rows(H) = count(i->any(!iszero,view(H,i,:)),axes(H,1))
function kernel_integer(A::AbstractMatrix{<:Integer})
    H,U = hnf(permutedims(A))
    permutedims(U[_nonzero_rows(H)+1:end,:])
end
function image_integer(A::AbstractMatrix{<:Integer})
    H = hnf(permutedims(A); transform=false)
    permutedims(H[1:_nonzero_rows(H),:])
end
"""Return nothing for an inconsistent integer system, otherwise x and a column kernel basis."""
function solve_integer(A::AbstractMatrix{<:Integer},b::AbstractVector{<:Integer})
    length(b)==size(A,1) || throw(DimensionMismatch("right-hand side"))
    H,U = hnf(permutedims(A))
    r = _nonzero_rows(H)
    residual = BigInt.(b)
    x = zeros(BigInt,size(A,2))
    for i in 1:r
        row = view(H,i,:)
        j = findfirst(!iszero,row)
        q,remnant = divrem(residual[j],row[j])
        iszero(remnant) || return nothing
        residual .-= q .* row
        x .+= q .* view(U,i,:)
    end
    all(iszero,residual) || return nothing
    (;x,kernel=permutedims(U[r+1:end,:]))
end
"""Characteristic polynomial of multiplication by h in the quotient by f."""
function charpoly(h::ModPoly,f::ModPoly)
    F=_check_poly_parents(h,f)
    n=degree(f)
    n>0 || throw(DomainError(f))
    g=rem(h,f)
    M=fill(zero(F),n,n)
    for j in 1:n
        for i in 1:n
            M[i,j]=coeff(g,i-1)
        end
        g=rem(LeftShift(g,1),f)
    end
    charpoly(M;ring=F)
end
charpoly(a::ExtensionElem)=charpoly(rep(a),modulus(a))
