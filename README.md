# NTL.jl

**A pure-Julia toolkit for exact finite-field arithmetic, number theory, polynomial computation, extension fields, exact linear algebra, generalized FFT, and lattice reduction.**

NTL.jl is inspired by C++ NTL, but it is not a wrapper or an API clone. It is designed to be approachable for new Julia users while providing efficient, tested implementations for research and teaching.

Version: **0.7.0**  
Tested with: **Julia 1.11.7**  
Compatibility target: **Julia 1.9+**

## Contents

- [Why NTL.jl?](#why-ntljl)
- [Installation](#installation)
- [A five-minute start](#a-five-minute-start)
- [Choosing the right field](#choosing-the-right-field)
- [Everyday operations](#everyday-operations)
- [Polynomials](#polynomials)
- [Multiplication and acceleration](#multiplication-and-acceleration)
- [GF2 and GF2X](#gf2-and-gf2x)
- [Extension fields](#extension-fields)
- [Zech fields for fast small-field arithmetic](#zech-fields-for-fast-small-field-arithmetic)
- [Matrices and exact linear algebra](#matrices-and-exact-linear-algebra)
- [LLL lattice reduction](#lll-lattice-reduction)
- [NTT and generalized FFT](#ntt-and-generalized-fft)
- [Number theory](#number-theory)
- [Important API reference](#important-api-reference)
- [Project architecture](#project-architecture)
- [Examples and tests](#examples-and-tests)
- [Comparison with C++ NTL](#comparison-with-c-ntl)
- [Comparison with GaloisFields.jl](#comparison-with-galoisfieldsjl)
- [Design principles](#design-principles)
- [Limitations](#limitations)
- [FAQ for beginners](#faq-for-beginners)
- [Contributing](#contributing)

## Why NTL.jl?

NTL.jl offers a coherent set of exact-mathematics tools rather than a loose collection of types.

- **Pure Julia:** no C++ compiler, NTL installation, binary artifact, or external Julia package is required for core functionality.
- **Exact arithmetic:** canonical residues, BigInt coefficients, exact polynomial operations, and verifiable lattice transformations.
- **Beginner-friendly parents:** `zzp{p}` stores a prime modulus in the type, while `ZZPField` and extension fields use explicit runtime parents.
- **Fast when appropriate:** machine-word multiplication, Karatsuba, generalized FFT, Zech discrete-log tables, and fast LLL with exact verification.
- **Safe fallbacks:** original schoolbook multiplication and exact LLL remain available. Fast paths never change mathematical results.
- **Tested invariants:** arithmetic is checked against independent BigInt calculations, polynomial identities, factor reconstruction, matrix identities, LLL conditions, and direct convolution.

The package includes more than finite fields. It also supports integer factorization, modular square roots, discrete logarithms, primitive roots, polynomial factorization, matrices, NTT, and LLL.

## Installation

### From a local checkout

```julia
using Pkg
Pkg.develop(path = "/path/to/NTLJulia")
using NTL
```

Replace `/path/to/NTLJulia` with the directory that contains `Project.toml`.

### Run the complete test suite

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

The current suite contains **4,335,143 passing assertions**. Runtime code depends only on Julia's `Random` standard library; `Test` is used only by the test target.

### HTML documentation

Open `docs/index.html` in a browser. The documentation is self-contained and can be read offline.

## A five-minute start

### 1. Create and use a prime field

```julia
using NTL

F = zzp(101)        # F is zzp{101}
a = F(120)          # 120 mod 101 = 19
b = F(-1)           # -1 mod 101 = 100

println(a + b)      # zzp{101}(18)
println(a * b)      # zzp{101}(24)
println(a / b)      # division in F_101
println(a^-1)       # inverse using gcdx and sign correction
println(a^1000)     # fast binary power

@assert rep(a) == 19
@assert Int(a) == 19
@assert modulus(a) == 101
@assert characteristic(F) == 101
@assert a * inv(a) == one(F)
@assert zero(F) == F(0)
@assert one(F) == F(1)
```

Every field element stores a canonical representative in `0:p-1`.

### 2. Iterate over a small field

```julia
for x in zzp{5}
    println(x)
end

@assert length(zzp{5}) == 5
@assert collect(zzp{5}) == zzp{5}.([0, 1, 2, 3, 4])
```

Elements are produced in canonical order, from `0` through `p-1`. Iteration is lazy; `collect(F)` allocates all `p` elements and should not be used for a very large field.

### 3. Multiply two polynomials

```julia
F = zzp(101)
x = xpoly(F)

f = x^3 + 2x + 1
g = x^2 + 3x + 2

@assert poly_mul(f, g) == f * g
@assert poly_mul(f, g; method = :schoolbook) == f * g
@assert poly_mul(f, g; method = :karatsuba) == f * g
@assert poly_mul(f, g; method = :fft) == f * g
```

Coefficients are stored **from low degree to high degree**. The zero polynomial has degree `-1`.

### 4. Solve a finite-field linear system

```julia
F = zzp(101)
A = F.([2 3 1; 4 1 5; 3 7 2])
b = F.([1, 2, 3])

x = solve(A, b)
@assert A * x == b

@assert det(A) != 0
@assert A * matrix_inverse(A) == F.([1 0 0; 0 1 0; 0 0 1])
```

### 5. Reduce an integer lattice basis

```julia
A = BigInt[1 1 1; -1 0 2; 3 5 6]
result = lll(A)  # fast by default, with exact verification

@assert result.transform * A == result.basis
@assert abs(det_exact(result.transform)) == 1
@assert result.verified
```

## Choosing the right field

| Need | Recommended parent | Why |
|---|---|---|
| Fast prime-field arithmetic with an `Int` modulus | `zzp(p)` | Type-level modulus, canonical `Int` representation, field iteration |
| Arbitrary-precision modulus, including a composite one | `ZZP(n)` | BigInt representation and explicit runtime parent |
| Binary field arithmetic | `GF2` | XOR/AND arithmetic and bit-packed polynomials |
| General polynomial extension | `extension_field(F, m)` | Polynomial quotient `F[x]/(m)` |
| Frequent arithmetic in a small extension | `zech_field(E)` | Discrete-log multiplication, inversion, and powers |

### `zzp` versus `ZZP`

```julia
F = zzp(101)             # machine-word prime field
R = ZZP(12)              # arbitrary-precision modular ring

@assert R(5) * R(5) == R(1)
@assert_throws DomainError inv(R(6))
```

`zzp` deliberately requires a prime modulus. `ZZP` also permits composite moduli, but only units are invertible.

## Everyday operations

### Construction and metadata

```julia
F = zzp(97)
a = F(-3)

rep(a)              # 94
Int(a)              # 94
modulus(a)          # 97
characteristic(a)   # 97
parent(a)           # zzp{97}
zero(a)             # zzp{97}(0)
one(a)              # zzp{97}(1)
iszero(a)           # false
isone(a)            # false
```

### Arithmetic

```julia
x + y
x - y
-x
x * y
x / y
inv(x)
x^e
```

Mixed integer arithmetic is reduced automatically:

```julia
@assert F(3) + 100 == F(0)
@assert 2F(3) == F(6)
@assert F(3) * big(8) == F(3)
```

### Random elements and vectors

```julia
using Random

rng = MersenneTwister(2026)
a = rand(rng, F)
v = rand(rng, F, 10)
w = VectorRandom(F, 4; rng = rng)
c = VectorConv(F, [-1, 0, 98, big(10)^40])
```

### Useful scalar helpers

| Function | Meaning |
|---|---|
| `rep(x)` | Canonical representative |
| `Int(x)` | Convert to machine integer |
| `BigInt(x)` | Convert to arbitrary-precision integer |
| `modulus(F)` | Modulus of the parent |
| `characteristic(F)` | Field characteristic |
| `parent(x)` | Parent field object |
| `zero(F)`, `one(F)` | Additive and multiplicative identities |
| `iszero(x)`, `isone(x)` | Identity tests |
| `add`, `sub`, `negate`, `mul`, `sqr`, `divide`, `power` | Procedural arithmetic aliases |
| `muladd(x,y,z)` | Fused `x*y+z` with fewer temporary objects |

## Polynomials

### Construction

```julia
F = zzp(101)
f = poly(F, [1, 2, 3])   # 1 + 2x + 3x^2
x = xpoly(F)              # x

@assert degree(f) == 2
@assert coeff(f, 0) == 1
@assert coeff(f, 2) == 3
@assert coeff(f, 20) == 0
@assert LeadCoeff(f) == 3
@assert ConstTerm(f) == 1
```

`poly(F, f)` is also valid for constructing from a compatible polynomial. `coefficients(f)` returns a copy of the coefficient vector.

### Core operations

| Function | Purpose |
|---|---|
| `DivRem(f,g)` | Polynomial quotient and remainder |
| `GCD(f,g)` / `gcd(f,g)` | Monic Euclidean GCD |
| `XGCD(f,g)` | Return `(d,s,t)` with `d = s*f + t*g` |
| `derivative(f)` / `diff(f)` | Formal derivative |
| `evaluate(f,x)` | Horner evaluation |
| `compose(f,g)` | Polynomial composition |
| `interpolate(F,xs,ys)` | Lagrange interpolation |
| `truncate(f,n)` | Keep terms of degree `< n` |
| `MulTrunc(f,g,n)` | Multiply then truncate |
| `InvTrunc(f,n)` | Newton truncated inverse |
| `MulMod(f,g,m)` | Product modulo polynomial `m` |
| `InvMod(f,m)` | Inverse in a quotient ring |
| `PowerMod(f,e,m)` | Binary powering modulo `m` |

### Factorization and irreducibility

```julia
F = zzp(97)
x = xpoly(F)
f = (x + 1)^3 * (x^2 + x + 1)

fac = factor(f)
@assert reconstruct(fac) == f
@assert all(is_irreducible(g) for (g, multiplicity) in fac.factors)
```

The factorization layer includes:

- square-free decomposition;
- distinct-degree factorization;
- Cantor–Zassenhaus equal-degree splitting;
- Rabin irreducibility testing;
- random irreducible polynomial generation with `BuildIrred`;
- root extraction with `FindRoots` and `FindRoot`.

## Multiplication and acceleration

Acceleration is **enabled by default**.

```julia
set_polynomial_acceleration!(true)
result = poly_mul(f, g)
```

The automatic strategy is:

1. Schoolbook for small degrees;
2. Karatsuba for medium degrees;
3. Generalized FFT for large degrees, using same-prime NTT when possible and compatible NTT primes with CRT otherwise.

### Explicit method selection

```julia
reference = poly_mul(f, g; method = :schoolbook)

@assert poly_mul(f, g; method = :karatsuba) == reference
@assert poly_mul(f, g; method = :ntt) == reference
@assert poly_mul(f, g; method = :fft) == reference
@assert poly_mul(f, g; method = :auto) == reference
```

To restore the original behavior everywhere:

```julia
set_polynomial_acceleration!(false)
```

### Why are there both `:ntt` and `:fft`?

- `:ntt` requests the same-prime transform and therefore requires the transform length to divide `p-1`.
- `:fft` is the generalized route. It can use compatible NTT primes and reconstruct exact integer coefficients with a coefficient bound or CRT.

If generalized FFT is unavailable for an input, automatic mode falls back to Karatsuba.

## GF2 and GF2X

`GF2` is the two-element field. `GF2X` stores a binary polynomial as bits of a `BigInt`, so bit `i` is the coefficient of `x^i`.

```julia
x = GF2X(2)             # x
f = x^2 + x + 1         # x^2 + x + 1

@assert is_irreducible(f)
@assert (x + 1)^2 == x^2 + 1  # characteristic two
@assert weight(f) == 3

fac = factor((x + 1)^3 * (x^2 + x + 1)^2)
@assert reconstruct(fac) == (x + 1)^3 * (x^2 + x + 1)^2
```

Useful functions include `bits`, `weight`, `NumBits`, `NumBytes`, `DivRem`, `XGCD`, `MulMod`, `PowerMod`, `factor`, `is_irreducible`, and `BuildIrredGF2X`.

## Extension fields

`extension_field(F, m)` constructs `F[x]/(m)`. The default constructor checks irreducibility.

```julia
F = zzp(5)
m = xpoly(F)^2 + 2
E = extension_field(F, m)
t = gen(E)

@assert t^2 == E(3)
@assert cardinality(E) == 25
@assert t * inv(t) == one(E)
@assert norm(t) == F(2)
@assert trace(t) == F(0)
```

For a binary example, `GF2E` is available:

```julia
x = GF2X(2)
aes_modulus = x^8 + x^4 + x^3 + x + 1
@assert is_irreducible(aes_modulus)

E = GF2E(aes_modulus)
a = gen(E) + 1

@assert cardinality(E) == 256
@assert a^255 == one(E)
@assert a * inv(a) == one(E)
```

`trace` and `norm` return values in the base field.

## Zech fields for fast small-field arithmetic

A Zech field is a discrete-log representation of the same extension field. It is not a different mathematical field.

For a field of order `q`, a primitive element `g` is chosen so that every nonzero element is `g^k`. Multiplication becomes addition of exponents modulo `q-1`; inversion and powering become sign change and scalar multiplication of exponents.

```julia
F = zzp(3)
E = extension_field(F, poly(F, [2, 1, 1]))
Z = zech_field(E)
x = ZechElem(Z, Z._generator)

@assert cardinality(Z) == 9
@assert x * inv(x) == one(Z)
@assert x^8 == one(Z)
@assert muladd(x, x, one(x)) == x * x + one(x)
@assert isbitstype(typeof(x))
```

### Why is it faster?

Ordinary extension multiplication performs polynomial multiplication and reduction. Zech multiplication performs table lookups:

- `log[a] + log[b]` for multiplication;
- `log[a] - log[b]` for division;
- `-log[a]` for inversion;
- `e * log[a]` for powers.

Use `extension_element(a)` to return to the ordinary polynomial extension representation.

### Limits and requirements

- The base must currently be a machine-word prime field `zzp`.
- The extension order must satisfy `q <= 65535`.
- A primitive multiplicative generator is required. `zech_field` searches for one automatically; an explicitly supplied non-generator is rejected.

## Matrices and exact linear algebra

The generic matrix routines work over supported fields.

```julia
F = zzp(101)
A = F.([2 3 1; 4 1 5; 3 7 2])
b = F.([1, 2, 3])

x = solve(A, b)
d = det(A)
invA = matrix_inverse(A)

@assert A * x == b
@assert A * invA == F.([1 0 0; 0 1 0; 0 0 1])

S = F.([1 2 3; 2 4 6; 0 1 1])
@assert iszero(det(S))
@assert rank(S) == 2
@assert all(iszero, S * nullspace(S)[1])
```

| Function | Purpose |
|---|---|
| `det(A)` | Determinant by Gaussian elimination |
| `det_exact(A)` | Exact integer determinant by Bareiss elimination |
| `rref(A)` | Reduced row-echelon form and pivot columns |
| `rank(A)` | Rank over an exact field |
| `solve(A,b)` | Solve a nonsingular square system |
| `matrix_inverse(A)` | Inverse over a field |
| `nullspace(A)` | Basis of the right nullspace |

## LLL lattice reduction

```julia
A = BigInt[1 1 1; -1 0 2; 3 5 6]
result = lll(A; mode = :fast, verify = true)

@assert result.transform * A == result.basis
@assert abs(det_exact(result.transform)) == 1
@assert result.verified
```

### Modes

| Mode | Behavior |
|---|---|
| `:fast` | Exact integer row operations, floating GSO for decisions, exact verification |
| `:exact` | Original exact-rational implementation |

The default fast mode verifies:

- `transform * A == basis`;
- `abs(det_exact(transform)) == 1`;
- exact size-reduction bounds;
- the exact Lovász condition.

If verification fails, the implementation automatically falls back to exact mode. `verify = false` skips the additional verification pass when that trade-off is intentional.

## NTT and generalized FFT

### Same-prime NTT

`ntt_convolve(a,b)` computes convolution over a supported prime field. A radix-two transform of length `n` requires `n` to be a power of two and to divide `p-1`.

```julia
F = zzp(97)  # 97-1 = 96, which is divisible by 32
a = rand(MersenneTwister(1), F, 32)
@assert intt(ntt(a)) == a

b = rand(MersenneTwister(2), F, 32)
@assert ntt_convolve(a, b) == convolve(a, b)
```

### Generalized FFT

`fft_convolve(a,b)` uses a same-prime fast path when possible. Otherwise it selects compatible NTT primes and reconstructs exact coefficients with a coefficient bound and, when needed, multi-prime CRT.

```julia
F = zzp(101)
f = random_poly(F, 1000)
g = random_poly(F, 1000)

reference = poly_mul(f, g; method = :schoolbook)
@assert poly_mul(f, g; method = :fft) == reference
```

The reconstruction is exact; it is not a floating-point approximation.

## Number theory

| Function | Purpose |
|---|---|
| `ProbPrime(n; rounds=16)` | Miller–Rabin primality test; deterministic in the signed `Int` range |
| `MillerWitness(n,w)` | Test whether `w` witnesses compositeness |
| `NextPrime(m)` | Least probable prime `>= m` |
| `PrevPrime(m)` | Greatest probable prime `<= m` |
| `GenPrime(bits)` | Generate a probable prime with an exact bit length |
| `RandomBits(bits)` | Random integer with at most `bits` bits |
| `RandomLen(bits)` | Random integer with exactly `bits` bits |
| `RandomBnd(n)` | Uniform random integer below `n` |
| `FactorInt(n)` | Factorization using Miller–Rabin and Pollard–Brent rho |
| `EulerPhi(n)` | Euler's totient from the factorization |
| `CRT(r,m)` | Chinese remaindering for pairwise-coprime moduli |
| `SqrRootMod(a,p)` | Tonelli–Shanks square root; returns `nothing` for a nonresidue |
| `DiscreteLog(a,b,p)` | Baby-step/giant-step logarithm |
| `PrimitiveRoot(p)` | Generator of the multiplicative group modulo a prime |

### A small example

```julia
@assert ProbPrime(97)
@assert !ProbPrime(341)
@assert NextPrime(100) == 101
@assert PrevPrime(100) == 97

fac = FactorInt(99991 * 10007)
@assert prod(p^e for (p, e) in fac; init = BigInt(1)) == 99991 * 10007

x, M = CRT([2, 3, 2], [3, 5, 7])
@assert (x, M) == (23, 105)

g = PrimitiveRoot(97)
y = powermod(g, 42, 97)
@assert powermod(g, DiscreteLog(g, y, 97), 97) == y
```

## Important API reference

### Prime-field and modular arithmetic

| Function | Purpose |
|---|---|
| `zzp(p)` | Validate a prime and return `zzp{p}` |
| `zzp_Min`, `zzp_Max` | Bounds of supported machine-word prime moduli |
| `AddMod`, `SubMod`, `NegateMod` | Overflow-safe modular add/subtract/negate on reduced `Int` values |
| `MulMod`, `SqrMod` | Modular multiplication and squaring |
| `InvMod`, `InvModStatus` | Inverse or a nonthrowing status result |
| `PowerMod` | Binary modular powering with negative-exponent support |
| `muladd(x,y,z)` | Fused `x*y+z` |

### Polynomial acceleration

| Function | Purpose |
|---|---|
| `poly_mul(f,g; method)` | Select schoolbook, Karatsuba, NTT, generalized FFT, or automatic mode |
| `set_polynomial_acceleration!(Bool)` | Enable/disable automatic acceleration |
| `polynomial_acceleration()` | Inspect the current mode |
| `fft_convolve(a,b)` | Generalized FFT with exact reconstruction |
| `ntt_convolve(a,b)` | Same-prime NTT convolution |
| `convolve(a,b)` | Direct schoolbook convolution |

### Extension fields and Zech fields

| Function | Purpose |
|---|---|
| `extension_field(F,m)` | Polynomial quotient extension |
| `GF2E(m)` | Binary extension |
| `gen(E)` | Polynomial indeterminate in the extension |
| `defining_polynomial(E)` | Defining polynomial |
| `extension_degree(E)` | Extension degree |
| `cardinality(E)` | Number of elements |
| `trace(a)`, `norm(a)` | Absolute trace and norm to the base field |
| `zech_field(E)` | Build a small-field discrete-log representation |
| `extension_element(a)` | Convert a Zech element back to an ordinary extension element |

### Matrix and lattice functions

| Function | Purpose |
|---|---|
| `det`, `det_exact`, `rank`, `rref`, `solve`, `matrix_inverse`, `nullspace` | Exact matrix algorithms |
| `gram_schmidt_exact(B)` | Exact rational Gram–Schmidt data |
| `lll(A; mode, verify)` | LLL with fast and exact modes |

## Project architecture

```text
NTLJulia/
├── Project.toml
├── README.md
├── docs/
│   └── index.html
├── src/
│   ├── NTL.jl              # Module entry, prime fields, low-level arithmetic
│   ├── zzp_big.jl          # Arbitrary-precision modular rings
│   ├── numbertheory.jl     # Primality, factorization, CRT, roots, discrete logs
│   ├── polynomial.jl       # Dense polynomials and factorization
│   ├── polynomial_mul.jl   # Schoolbook, Karatsuba, and acceleration strategies
│   ├── binary.jl           # GF2 and bit-packed GF2X
│   ├── extension.jl        # Polynomial extension fields
│   ├── zech.jl             # Isbits Zech representation and tables
│   ├── matrix.jl           # Determinants, RREF, solving, rank, nullspace
│   ├── lattice.jl          # Original exact LLL implementation
│   ├── lll_fast.jl         # Fast LLL with exact verification and fallback
│   ├── transform.jl        # Same-prime NTT and optimized butterflies
│   └── fft.jl              # Generalized FFT and CRT reconstruction
├── examples/
│   ├── 01_basic_arithmetic.jl
│   ├── 02_vectors_and_random.jl
│   ├── 03_linear_system.jl
│   ├── 04_number_theory.jl
│   ├── 05_polynomial_factorization.jl
│   ├── 06_extension_gf2.jl
│   ├── 07_matrix_lattice.jl
│   ├── 08_ntt_convolution.jl
│   ├── 09_polynomial_acceleration.jl
│   ├── 10_generalized_fft.jl
│   ├── 11_fast_lll.jl
│   └── 12_zech_extension.jl
└── test/
    ├── runtests.jl
    ├── zzp_numtheory.jl
    ├── polynomial.jl
    ├── binary_extension.jl
    ├── matrix_lattice_transform.jl
    ├── acceleration.jl
    ├── lll_fast.jl
    └── zech.jl
```

The public API is grouped by mathematical role rather than by C++ header names, so the implementation is easier to navigate for new users.

## Examples and tests

Run any example with:

```bash
julia --project=. examples/01_basic_arithmetic.jl
```

| Script | Demonstrates |
|---|---|
| `01_basic_arithmetic.jl` | Prime-field construction, arithmetic, low-level helpers, iteration |
| `02_vectors_and_random.jl` | Random elements, conversion, vector operations |
| `03_linear_system.jl` | Solving a finite-field linear system |
| `04_number_theory.jl` | Primality, factorization, CRT, modular roots, discrete logs, prime generation |
| `05_polynomial_factorization.jl` | Factor reconstruction, irreducibility, interpolation, modular powers |
| `06_extension_gf2.jl` | GF(2^8), prime-field extensions, extension-field polynomials |
| `07_matrix_lattice.jl` | Matrix solving, inverse, nullspace, exact LLL |
| `08_ntt_convolution.jl` | NTT round trips and convolution |
| `09_polynomial_acceleration.jl` | Schoolbook, Karatsuba, NTT, and automatic selection |
| `10_generalized_fft.jl` | Generalized FFT and exact coefficient reconstruction |
| `11_fast_lll.jl` | Fast and exact LLL modes with verification |
| `12_zech_extension.jl` | Zech elements, inverse, powers, and `muladd` |

## Comparison with C++ NTL

The following measurements are **workload-specific examples**, not universal performance guarantees. They were collected on a 2-core server using official C++ NTL 11.4.2 with native SIMD and Julia 1.11.7.

| Operation | C++ NTL (ms) | NTL.jl (ms) | Julia / C++ |
|---|---:|---:|---:|
| Small modular multiplication | 8.5685 | 20.0726 | 2.34x |
| BigInt modular exponentiation | 0.0101 | 0.0189 | 1.87x |
| Polynomial multiplication | 0.1025 | 2.5212 | 24.60x |
| NTT convolution | 0.8374 | 5.8378 | 6.97x |
| Polynomial factorization | 0.0297 | 2.0918 | 70.43x |
| Finite-field determinant | 0.0395 | 0.4852 | 12.28x |
| LLL | 0.0720 | 2.0970 | 29.13x |
| Primality test | 0.0680 | 0.0295 | 0.43x |

### What these results mean

- Scalar modular arithmetic is already practical for many workloads.
- C++ NTL remains faster for large FFT-based multiplication, factorization, and LLL because it has specialized tables, precomputation, pointer-level loops, and long-standing native optimizations.
- Julia's generalized FFT substantially reduced the polynomial multiplication gap compared with the original schoolbook implementation.
- Fast LLL reduced the LLL-specific slowdown from over 1,500x to about 29x while retaining exact verification.
- The primality result is not evidence that Julia is generally faster: the compared witness sets and rounds differ.

## Comparison with GaloisFields.jl

GaloisFields.jl is particularly strong for small finite fields and extension arithmetic. It uses Zech logarithms for fast multiplication and has a convenient macro interface.

| Operation | NTL.jl (ms) | GaloisFields (ms) |
|---|---:|---:|
| Prime-field mul/add x1e6, p=101 | 48.15 | 53.54 |
| Prime-field fused muladd x1e6 | 26.24 | 53.54 |
| Prime-field power x10000, p=101 | 0.718 | 0.687 |
| Prime-field inverse x10000, p=101 | 0.570 | 0.577 |
| Zech extension ops x10000 | 1.283 | 1.191 |
| Zech muladd x10000 | 1.064 | 1.191 |
| Zech inverse x10000 | 0.660 | 0.746 |
| Zech power x10000 | 0.783 | 0.680 |
| Polynomial multiplication, deg 999 | 1.206 | 6.347 |

### Interpretation

- Prime-field arithmetic is comparable, and fused `muladd` is faster in the tested loop.
- Zech arithmetic is also comparable. The small differences are implementation overheads rather than fundamental algorithmic differences.
- NTL.jl is faster for the tested large polynomial multiplication because of generalized FFT and exact coefficient reconstruction.
- GaloisFields remains a good choice for compact small-field code and extension-tower syntax; NTL.jl covers a broader set of algorithms including matrices, LLL, generalized FFT, and number theory.

## Design principles

### Mathematical correctness before speed

Fast paths are checked against the original algorithms. A user can always request `:schoolbook` or `mode = :exact` to obtain the original behavior.

### Explicit mathematical parents

`zzp{p}` keeps the modulus in the type, so different moduli cannot be mixed accidentally. `ZZPField`, `ExtensionField`, and `ZechField` make runtime parent relationships explicit.

### Canonical and exact representations

Prime-field elements are always stored in `0:p-1`. Polynomial coefficients are stored low degree first. Zech elements use isbits integer representatives with a parent-table identifier.

### Acceleration with safe fallbacks

The package prefers degree-based and representation-based acceleration rather than one universal algorithm. Unsupported NTT lengths, large coefficient bounds, or failed LLL verification trigger a safe fallback.

### Verification as part of the implementation

The LLL fast path checks the returned basis and transformation with exact rational GSO. Polynomial acceleration is tested by comparing every result with schoolbook and direct convolution.

## Limitations

- This is not a full replacement for every NTL module. Arbitrary-precision real numbers (`RR`), HNF, BKZ, complete integer-polynomial factorization, and C++ ABI compatibility are not included.
- Radix-two NTT requires a power-of-two transform length and an appropriate root of unity.
- Generalized FFT currently supports scalar modular coefficients; extension-field coefficients require other multiplication paths.
- Zech tables require a primitive generator and are limited to fields with at most 65,535 elements.
- Fast LLL is intended for reproducible exact arithmetic on modest bases, not as a replacement for NTL's highly tuned large-basis reduction routines.
- Polynomial factorization and matrix algorithms are not claimed to match C++ NTL's performance for very large inputs.

## FAQ for beginners

### Why does `zzp(15)` fail?

`zzp` represents a prime field. Fifteen is composite. Use `ZZP(15)` if you intentionally need a modular ring, and remember that only units have inverses.

### Why is `F(2)` zero when `F = zzp(2)`?

Construction always reduces modulo the field modulus: `2 mod 2 = 0`.

### Why can I not add `zzp{7}(1)` to `zzp{11}(1)`?

They belong to different fields. Choose one parent explicitly instead of relying on an implicit and potentially incorrect conversion.

### Why does `collect(F)` use memory even though iteration is lazy?

Iteration produces elements one at a time, but `collect` allocates a vector of all `p` elements. Use `length(F)` for size and `for x in F` for lazy traversal.

### Why does NTT report that `n` must divide `p-1`?

A primitive `n`-th root of unity exists in `F_p` only when `n` divides `p-1`. If that condition fails, automatic multiplication uses generalized FFT or Karatsuba instead.

### Why does generalized FFT sometimes use other primes?

It allows convolution when the original field lacks a suitable root. Exact integer bounds and CRT ensure that the reconstructed result is correct before reduction to the original modulus.

### When should I use `zech_field` instead of `extension_field`?

Use it for frequent multiplication, inversion, and powers in a small extension. It builds a discrete-log table once and then performs fast table-based operations. For large extensions or polynomial-heavy theory, use the ordinary extension representation.

### Is `verify = false` in LLL unsafe?

It does not change the integer arithmetic of the basis or transformation, but it skips the additional exact verification step. Use it only when you intentionally want to avoid that cost.

## Contributing

Contributions should preserve the package's verification-first style.

1. Add mathematical identities and regression tests before optimizing.
2. Compare with schoolbook, direct convolution, or C++ NTL where practical.
3. Keep the original algorithm available as an explicit fallback.
4. Document supported moduli, parent requirements, complexity, and failure modes.
5. Run the full test suite and all examples before pushing.
6. Update `README.md`, `docs/index.html`, and relevant examples for every public interface change.

Useful invariants include `U*A == B`, `abs(det_exact(U)) == 1`, factor reconstruction, `A*x == b`, and exact coefficient reconstruction.
