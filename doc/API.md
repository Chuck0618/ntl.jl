+ coeff(zz_pX{T}, i ::Int)::zz_p{T}
+ coeff(zz_pX{T})

Get the coefficient of x^i. 

---
+ normalize!(zz_pX{T})

Delete the zero leading terms.

---

+ iszero(zz_pX)::bool
+ iszero(zz_p)::bool

Check whether the value is zero

---
+ zz_pX(T::Int)::DataType

Construct the type: polynomial ring over zz_p;

---