export bkz
# Enumerate exactly the integer combinations in the projected block whose
# squared length is strictly smaller than that of its first vector.
function _shorter_block(mu,norms,lo,hi)
    z = zeros(BigInt,hi-lo+1)
    bound = norms[lo]
    answer = nothing
    function visit(i,used)
        if i<lo
            if used>0
                bound = used
                answer = copy(z)
            end
            return
        end
        center = -sum((mu[j,i]*z[j-lo+1] for j in i+1:hi);init=big(0)//1)
        radius2 = (bound-used)/norms[i]
        radius2>0 || return
        radius = isqrt(fld(numerator(radius2),denominator(radius2)))+1
        for a in ceil(BigInt,center)-radius:floor(BigInt,center)+radius
            cost = used+(a-center)^2*norms[i]
            cost<bound || continue
            z[i-lo+1]=a
            visit(i-1,cost)
        end
    end
    visit(hi,big(0)//1)
    answer
end

"""Exact projected-block SVP enumeration; exponential cost in block_size."""
function bkz(A::AbstractMatrix{<:Integer};block_size::Integer=20,delta=99//100,
             verify::Bool=true,maxiter::Integer=100000)
    n,m=size(A)
    1<=n<=m && block_size>=2 && maxiter>0 || throw(ArgumentError("expected a nonempty independent row basis, block_size ≥ 2 and positive maxiter"))
    d=_lll_delta(delta)
    B,U=_lll_exact_legacy(A;delta=d,maxiter=maxiter)
    for _ in 1:maxiter
        mu,norms=gram_schmidt_exact(B)
        changed=false
        for lo in 1:n-1
            hi=lo+min(block_size,n-lo+1)-1
            z=_shorter_block(mu,norms,lo,hi)
            isnothing(z) && continue
            # V*z=e₁; transpose(inv(V)) completes z to a unimodular row basis.
            H,V=hnf(reshape(z,:,1))
            H[1,1]==1 || error("shortest coefficient vector must be primitive")
            T=permutedims(BigInt.(matrix_inverse(Rational{BigInt}.(V))))
            I=lo:hi
            B[I,:]=T*B[I,:]
            U[I,:]=T*U[I,:]
            B,W=_lll_exact_legacy(B;delta=d,maxiter=maxiter)
            U=W*U
            changed=true
            break
        end
        if !changed
            verify && !_lll_verify(A,B,U,d) && error("invalid exact LLL certificate")
            return (;basis=B,transform=U,verified=verify)
        end
    end
    error("BKZ tour limit exceeded")
end
