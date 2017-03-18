DoubleDoubles.jl
===============

(The module was renamed `DoubleDoubles` as it provides the type `DoubleDouble`.)

This is work in progress. Don't switch from DoubleDouble.jl just yet. A lot of the code is completely untested.

`DoubleDoubles.jl` is a Julia package for performing extended-precision arithmetic using pairs of floating-point numbers. This is commonly known as "double-double" arithmetic, as the most common format is a pair of C-doubles (`Float64` in Julia), although `Double`s will actually work for any floating-point type, including itself. Its aim is to provide accurate results without the overhead of `BigFloat` types.

The core routines are based on the ideas and algorithms of [Dekker (1971)][dekker1971].

Interface
---------
The main type is `Double`, with two floating-point fields: `hi`, storing the leading bits, and `lo` storing the remainder. `hi` is stored to full precision and rounded to nearest.  When creating a `Double` directly using the type constructor, the user is responsible for ensuring that `abs(x.lo) <= 0.5 * eps(x.hi)`.
However, `Double`s are typically created by conversion from other numeric types.

```julia
julia> using DoubleDoubles

julia> x = DoubleDouble(pi)
3.14159265358979323846264338327953

julia> (x.hi, x.lo)
(3.141592653589793,1.2246467991473532e-16)

julia> eps(x.hi)
4.440892098500626e-16

julia> xx = QuadDouble(pi)
3.1415926535897932384626433832795028841971693993751058209749445924

julia> (xx.hi.hi, xx.hi.lo, xx.lo.hi, xx.lo.lo)
(3.141592653589793,1.2246467991473532e-16,-2.764667309787496e-32,1.4797009536535408e-48)
```

The other type defined is `Single`, which is a `Double` where the `lo` field is guaranteed to be zero.
(Attempting a conversion to `Single` will throw an `InexactError` if this is not the case.)
Operations on `Single` typically return `Double` and will often be faster than the corresponding operations on `Double`.


### Note
Constructs like `Double(pi)` that used to default to `Double{Float64}(pi)` are depreciated.
Either write out `Double{Float64}` explicitly, or use the type alias `DoubleDouble`.

`Double{BigFloat}` is not allowed. The base type must be immutable. (Use `setprecision` to get higher-precision `BigFloat`s.)

Examples
---------
### Exact products and remainders

By exploiting this property, we can compute exact products of floating point numbers.

```julia
julia> u, v = 64 * rand(), 64 * rand()
(44.438125149181445, 46.96067529434315)

julia> w = Single(u) * Single(v)
2.08684436582009398756273247966684e+03

julia> (w.hi, w.lo)
(2086.844365820094, -3.871317795744025e-14)
```
Note that the product of two `Single`s is a `Double`: the `hi` element of this
double is equal to the usual rounded product, and the `lo` element contains the exact
difference between the exact value and the rounded.

This can be used to get an accurate remainder
```julia
julia> r = rem(w, 1.0)
8.44365820093987562732479666836822e-01

julia> Float64(r)
0.8443658200939875
```

This is much more accurate than taking ordinary products, and gives the same answer as using `BigFloat`s:
```julia
julia> rem(u*v, 1.0)
0.8443658200940263

julia> Float64(rem(big(u) * big(v), 1.0))
0.8443658200939875
```
However, since the `DoubleDouble` version is carried out using ordinary floating-point operations, it is of the order of 1000x faster than the `BigFloat` version.

### Correct rounding with non-exact floats

If a number cannot be exactly represented by a floating-point number, it may be rounded incorrectly when used later, e.g.
```julia
julia> pi * 0.1
0.3141592653589793

julia> Float64(big(pi) * 0.1)
0.31415926535897937
```
We can also do this computation using `Double`s (note that the promotion rules mean that only one needs to be specified):
```julia
julia> Float64(DoubleDouble(pi) * 0.1)
0.31415926535897937

julia> Float64(pi * Single(0.1))
0.31415926535897937
```

### Emulated FMA

The [fused multiply-add (FMA)](http://en.wikipedia.org/wiki/Multiply%E2%80%93accumulate_operation) operation is an intrinsic floating-point
operation that allows the evaluation of `a * b + c`, with rounding occurring only
at the last step. This operation is unavailable on 32-bit x86 architecture, and available
only on the most recent x86_64 chips, but can be emulated via double-double arithmetic:

```julia
julia> 0.1 * 0.1 + 0.1
0.11000000000000001

julia> Float64(big(0.1) * 0.1 + 0.1)
0.11

julia> Base.fma(a::Float64,b::Float64,c::Float64) = Float64(Single(a) * Single(b) + Single(c))
fma (generic function with 1 method)

julia> fma(0.1, 0.1, 0.1)
0.11
```

[dekker1971]: http://link.springer.com/article/10.1007%2FBF01397083  "T.J. Dekker (1971) 'A floating-point technique for extending the available precision', Numerische Mathematik, Volume 18, Issue 3, pp 224-242"
