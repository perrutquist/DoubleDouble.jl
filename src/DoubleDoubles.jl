module DoubleDoubles

using Zeros

export AbstractDouble, Double, DoubleDouble, QuadDouble, OctaDouble,
       Single, SingleSingle, SingleSingleSingle

import Base:
    convert,
    *, +, -, /, sqrt, <, >, <=, >=, ==, !=, ldexp,
    rem, rand, promote_rule, copysign, flipsign, abs,
    round, floor, ceil, trunc, div, maxintfloat, typemin, typemax,
    show, big, precision, realmin, realmax, exponent, eps

abstract type AbstractDouble{T} <: AbstractFloat end

"A Single is a Double where the lower part is exactly zero."
struct Single{T<:AbstractFloat} <: AbstractDouble{T}
    hi::T
    lo::Zero
end
Single(hi::AbstractFloat) = Single(hi,Zero())

#"Just like error(), but not in-lined."
#@noinline throwerror(s...) = throw(ErrorException(s...))

#"Test if a type is immutable"
#@noinline isimmutable_type{T}(::Type{T}) = isimmutable(zero(T))

@generated checkimmutable(hi) = hi.mutable ? :( throw("Needed an immutable type!") ) : nothing

"Double{T} uses two instances of type T to represent a floating point number."
# In a Double, hi uses the full mantissa, and abs(lo) <= 0.5eps(hi)
# This must be enforced by the method that creates the double.
# No check at runtime for performance reasons.
struct Double{T<:AbstractFloat} <: AbstractDouble{T}
    hi::T
    lo::T
    function Double{T}(hi::T,lo::T) where T<:AbstractFloat
      checkimmutable(hi)
      if T<:BigFloat
        error("The parts of a Double must be immutable. BigFloat is not allowed")
      end
      if T<:Single # Tests at compile-time,
        error("Creating a Double{Single{...}} makes no sense.")
      end
      if T==Float16 # Tests at compile-time,
        error("Double{Flat16} not supported at the moment. (Use Float32 instead.)")
      end
      # TODO: Test if precision is too high, potentially causing underflow in .lo
      new(hi,lo)
    end
end

# Allow conversion via e.g. Double{Float64}(x)
Double{T}(x::Real) where T<:AbstractFloat = convert(Double{T}, x)

#Is this not redundant?
Double(x::T,y::T) where T<:AbstractFloat = Double{T}(x,y)

function Double(x::Real)
  warn("Double(x) is depreciated. Use Double{Float64}(x) or DoubleDouble(x) to convert x to Double{Float64}.")
  Double{Float64}(x)
end

# A Double with a Zero as lo is by definition a Single.
Double(hi::AbstractFloat, ::Zero) = Single(hi)
Double{T}(hi::AbstractFloat, ::Zero) where T = Single{T}(convert(T, hi))

"Convert Single to Double"
nosingle(x) = x
nosingle(x::Double) = Double(nosingle(x.hi), nosingle(x.lo))
nosingle(x::Single) = Double(nosingle(x.hi), nosingle(zero(x.hi)))

"Translate Single type into the corresponding Double type"
nosingle(T::Type) = T
# We don't bother with Double, because Double{Single} is not allowed.
nosingle(::Type{Single{T}}) where {T<:AbstractFloat} = Double{nosingle(T)}

# "Normalize" Doubles to ensure abs(lo) <= 0.5eps(hi)
# assumes abs(u) > abs(v): if not, use Single + Single
function normalize_double(u::T, v::T) where T<:AbstractFloat
  w = u + v
  Double{T}(w, (u-w) + v)
end
normalize_double(u::T, v::T) where T<:Single = normalize_double(nosingle(u), nosingle(v))
normalize_double(u::AbstractFloat, ::Zero) = Single(u)

# This is a lower bound on the precision. If the binary representation of a
# number has a lot of zeros in the middle, then higher precision is achieved.
# For example, 1 + 2^-300 can be represented exactly by a Double{Float64}
precision(::Type{Double{T}}) where T<:AbstractFloat = 2*precision(T)

# The precision of Single is defined the same as the precision of the
# correspoinding Double. The digits "stored" in the always-zero .lo
# field count! (This is why Single contructors may throw InexactError.)
precision(::Type{Single{T}}) where T<:AbstractFloat = 2*precision(T)

# ldexp can be done without renormalization
ldexp(x::AbstractDouble, n::Integer) = Double{typeof(x.hi)}(ldexp(x.hi,n), ldexp(x.lo,n))

"""
    ldmul(x,poweroftwo)

Returns x*poweroftwo, in the type of x.
May give incorrect results if poweroftwo is not a power of two.
For x::Double, this is faster than normal multiplication because no
renormalization is needed.
"""
ldmul(x::AbstractFloat, poweroftwo::Real) = convert(typeof(x), x*poweroftwo)
function ldmul(x::AbstractDouble, poweroftwo::Real)
  Double(ldmul(x.hi,poweroftwo),ldmul(x.lo,poweroftwo))
end

"""
A number such that half the significand bits are used for the integer part
and the other half of the bits are used for the fractional part.
"""
@generated halfprec_constant(::Type{T}) where T<:AbstractFloat = convert(T,exp2(ceil(precision(T)/2)))

"round floats to half-precision"
# TODO: fix overflow for large values
halfprec(x) = (p = x*halfprec_constant(typeof(x)); (x-p)+p)

"Split a float into a two-tuple of its lower and upple half"
function splitprec(x::AbstractFloat)
    h = halfprec(x)
    h, x-h
end

function splitprec(x::Double)
    Single(x.lo), Single(x.hi)
end

## Conversion to Single
convert(::Type{Single{T}}, x::T) where T<:AbstractFloat = Single(x)
convert(::Type{Single{T}}, x::Double{T}) where T<:AbstractFloat = x.lo == zero(T) ? Single(x.hi) : throw(InexactError())
convert(::Type{Single{T}}, x::Single{T}) where T<:AbstractFloat = x # needed because Single <: AbstractFloat

function singleconvert(::Type{T}, x::Real) where T<:AbstractFloat
  z = convert(T,x)
  if convert(typeof(x),z) != x #debug
      println("difference: ", repr(x-convert(typeof(x),z)))
      println("original: ", repr(x))
      println("converted: ", repr(z))
      println("converted back: ", repr(convert(typeof(x),z)))
  end
  convert(typeof(x),z) == x || throw(InexactError(:singleconvert, T, x))
  Single(z)
end
convert(::Type{Single{T}}, x::Real) where T<:AbstractFloat = singleconvert(T, x)
# Avoid ambiguity with Rational-to-AbstractFloat conversion in rational.jl
convert(::Type{Single{T}}, x::Rational) where T<:AbstractFloat = singleconvert(T, x)

## Conversion to Double
function doubleconvert(::Type{T}, x::Real) where T<:AbstractFloat
  z = convert(T, x)
  y = convert(typeof(x), z)
  Double{T}(z, convert(T, x-y))
end
function doubleconvert(::Type{T}, x::Unsigned) where T<:AbstractFloat
  z = convert(T, x)
  y = convert(typeof(x), z)
  if y>x
    return Double{T}(z, -convert(T, y-x))
  else
    return Double{T}(z, convert(T, x-y))
  end
end

convert(::Type{Double{T}}, x::T) where T<:AbstractFloat = Double{T}(x, zero(T))
convert(::Type{Double{T}}, x::Single{T}) where T<:AbstractFloat = Double{T}(x.hi, zero(T))
convert(::Type{Double{T}}, x::Rational) where T<:AbstractFloat = Double{T}(x.num)/Double{T}(x.den)
convert(::Type{Double{T}}, x::Double{T}) where T<:AbstractFloat = x # needed because Double <: AbstractFloat

# Disambiguation
convert(::Type{Double{T1}}, x::AbstractDouble{T2}) where T2<:AbstractFloat where T1<:AbstractFloat =
  doubleconvert(T1, x)
convert(::Type{Single{T}}, x::Single) where T<:AbstractFloat = singleconvert(T, x)
convert(::Type{Single{T}}, x::Double) where T<:AbstractFloat = singleconvert(T, x)
convert(::Type{Double{T}}, x::Real) where T<:AbstractFloat = doubleconvert(T, x)
convert(::Type{Double{T}}, x::Zero) where T<:AbstractFloat = Double{T}(0)
convert(::Type{Single{T}}, x::Zero) where T<:AbstractFloat = Single{T}(0)

function convert(::Type{BigFloat}, x::AbstractDouble{T}) where T<:AbstractFloat
  setprecision(BigFloat, 2*precision(T)) do
    convert(BigFloat, x.hi) + convert(BigFloat, x.lo)
  end
end

# Macro will do the BigFloat generation of irrationals at compile-time.
@generated function doublesym(::Type{Double{T}}, ::Irrational{sym}) where {T<:AbstractFloat,sym}
  Double{T}(setprecision(()->BigFloat(eval(sym)),BigFloat,2*precision(T)+32))
end
convert(::Type{Double{T}}, sym::Irrational) where T<:AbstractFloat = doublesym(Double{T}, sym)

## Conversion from Single/Double
convert(::Type{T}, x::AbstractDouble{T}) where T<:AbstractFloat = x.hi
convert(::Type{T1}, x::AbstractDouble{T2}) where T2<:AbstractFloat where T1<:AbstractFloat =
convert(T1, x.hi) + (precision(T1) > precision(T2) ? convert(T1, x.lo) : 0)

function convert(::Type{T1}, x::AbstractDouble{T2}) where T2<:AbstractFloat where T1<:Integer
  # TODO: Check range
  #if x.hi > convert(baseof(typeof(x)), typemax(T1)) || x.hi < convert(baseof(typeof(x)), typemin(T1))
  #  throw(InexactError)
  #end
  if T1<:Unsigned && x.lo<0
    return convert(T1, x.hi) - convert(T1, -x.lo)
  else
    return convert(T1, x.hi) + convert(T1, x.lo)
  end
end

## Disambiguation in convert Double to Double
function convert(::Type{Double{T}}, x::Double) where T<:AbstractFloat
  if typeof(x) == T
    Double{T}(x,zero(T))
  else
    if precision(x.hi) > 2*precision(T)
      Double{T}(x.hi)
    else
      Double{T}(x.hi) + Double{T}(x.lo) # Faster with normalize_double instead of plus ?
    end
  end
end

## promotion
promote_rule(::Type{Double{T}}, ::Type{T2}) where T2<:Integer where T<:AbstractFloat = Double{T}
promote_rule(::Type{Single{T}}, ::Type{T2}) where T2<:Integer where T<:AbstractFloat = Double{nosingle(T)}

promote_rule(::Type{Single{T1}}, ::Type{T2}) where T2<:AbstractFloat where T1<:AbstractFloat = 2*precision(T1) > precision(T2) ? Double{nosingle(T1)} : nosingle(T2)
promote_rule(::Type{Double{T1}}, ::Type{T2}) where T2<:AbstractFloat where T1<:AbstractFloat = 2*precision(T1) > precision(T2) ? Double{T1} : nosingle(T2)
promote_rule(::Type{Double{T}}, ::Type{Single{T}}) where T<:AbstractFloat = Double{T}

promote_rule(::Type{Irrational{s}}, ::Type{Single{T}}) where T<:AbstractFloat where s<:Integer = Double{nosingle(T)}
promote_rule(::Type{Irrational{s}}, ::Type{Double{T}}) where T<:AbstractFloat where s<:Integer = Double{T}

promote_rule(::Type{AbstractDouble{T}}, ::Type{BigFloat}) where T<:AbstractFloat = BigFloat

# Type aliases

Quad{T} = Double{Double{T}}
Octa{T} = Double{Quad{T}}

const DoubleDouble = Double{Float64}
const QuadDouble = Quad{Float64}
const OctaDouble = Octa{Float64}

SingleSingle{T} = Single{Single{T}} # A Quad with .hi.hi as only nonzero
SingleSingleSingle{T} = Single{SingleSingle{T}} # An Octa

# Comparisons

for op in (:<, :>)
  @eval $op(x::AbstractDouble{T},y::AbstractDouble{T}) where T<:AbstractFloat =
    $op(x.hi,y.hi) || (x.hi == y.hi && $op(x.lo,y.lo))
end
for (op, nop) in ((:<=, :>), (:>=, :<))
  @eval $op(x::AbstractDouble{T},y::AbstractDouble{T}) where T<:AbstractFloat =
    !$nop(x,y)
end

==(x::AbstractDouble{T},y::AbstractDouble{T}) where T<:AbstractFloat = (x.hi==y.hi && x.lo==y.lo)
!=(x::AbstractDouble{T},y::AbstractDouble{T}) where T<:AbstractFloat = (x.hi!=y.hi || x.lo!=y.lo)
# hash(x::Single, h::UInt) = hash(hash(x.hi),h);
# hash(x::Double, h::UInt) = hash(hash(x.hi),hash(x.lo,h));

# TODO eliminate branches
# add12
function +(x::Single{T},y::Single{T}) where T<:AbstractFloat
    abs(x.hi) > abs(y.hi) ? normalize_double(x.hi, y.hi) : normalize_double(y.hi, x.hi)
end

# Dekker add2
function +(x::AbstractDouble{T}, y::AbstractDouble{T}) where T<:AbstractFloat
    r = x.hi + y.hi
    s = abs(x.hi) > abs(y.hi) ? (((x.hi - r) + y.hi) + y.lo) + x.lo : (((y.hi - r) + x.hi) + x.lo) + y.lo
    normalize_double(r, s)
end

-(x::AbstractDouble) = Double(-x.hi, -x.lo)

function -(x::Single{T}, y::Single{T}) where T<:AbstractFloat
    abs(x.hi) > abs(y.hi) ? normalize_double(x.hi, -y.hi) : normalize_double(-y.hi, x.hi)
end

function -(x::AbstractDouble{T}, y::AbstractDouble{T}) where T<:AbstractFloat
    r = x.hi - y.hi
    s = abs(x.hi) > abs(y.hi) ? (((x.hi - r) - y.hi) - y.lo) + x.lo : (((-y.hi - r) + x.hi) + x.lo) - y.lo
    normalize_double(r, s)
end

# TODO FMA version
# Dekker mul12
function *(x::Single{T},y::Single{T}) where T<:AbstractFloat
    hx,lx = splitprec(x.hi)
    hy,ly = splitprec(y.hi)
    z = x.hi*y.hi
    normalize_double(z, ((hx*hy-z) + hx*ly + lx*hy) + lx*ly)
end

# Dekker mul2
function *(x::AbstractDouble{T}, y::AbstractDouble{T}) where T<:AbstractFloat
    c = Single(x.hi) * Single(y.hi)
    cc = (x.hi * y.lo + x.lo * y.hi) + c.lo
    normalize_double(c.hi, cc)
end

# Dekker div2
function /(x::AbstractDouble{T}, y::AbstractDouble{T}) where T<:AbstractFloat
    c = x.hi / y.hi
    u = Single(c) * Single(y.hi)
    cc = ((((x.hi - u.hi) - u.lo) + x.lo) - c*y.lo)/y.hi
    normalize_double(c, cc)
end

# Dekker sqrt2
function sqrt(x::AbstractDouble{T}) where T<:AbstractFloat
    if x.hi <= 0
        throw(DomainError("sqrt will only return a complex result if called with a complex argument."))
    end
    c = sqrt(x.hi)
    u = Single(c)*Single(c)
    cc = (((x.hi - u.hi) - u.lo) + x.lo)*convert(T,0.5)/c
    Double{T}(c, cc)
end

"The highest underlying primitive float ( x.hi.hi.hi... )"
highfloat(x::AbstractFloat) = x
highfloat(x::AbstractDouble) = highfloat(x.hi)

"High accuracy remainder."
rem_accurate(x::T,d::T) where {T<:AbstractFloat} = begin
    x = Single(x)
    d = Single(d)
    T(x-d*div(x,d))
end
rem(x::AbstractDouble{T},d::AbstractDouble{T}) where {T<:AbstractFloat} = x-d*div(x,d)
rem(x::AbstractDouble,d::Real) = normalize_double(rem(x.hi,d), rem(x.lo,d))

signbit(x::AbstractDouble) = signbit(x.hi)

for op in [:round, :floor, :ceil, :trunc]
  @eval $op(x::AbstractDouble) = abs(x.lo) >= 1 ? Double(x.hi, $op(x.lo)) : Double($op(x.hi),zero(x.lo))
end
div(x::AbstractDouble{T},y::AbstractDouble{T}) where T<:AbstractFloat = trunc(x/y)

# TODO: Test if the below is really better than the default impelentation
for op in [:copysign, :flipsign]
  @eval $op(x::AbstractFloat,y::AbstractDouble) = $op(x,highfloat(y))
  @eval $op(x::AbstractDouble,y::AbstractFloat) = Double($op(x.hi,y), $op(x.lo,y))
  @eval $op(x::AbstractDouble,y::AbstractDouble) = Double($op(x.hi,highfloat(y)), $op(x.lo,highfloat(y)))
end
abs(x::AbstractDouble) = flipsign(x,x)

realmin(::Type{Single{T}}) where T<:AbstractFloat = Single{T}(ldexp(realmin(T),precision(T)))
realmin(::Type{Double{T}}) where T<:AbstractFloat = Double{T}(ldexp(realmin(T),precision(T)))

exponent(x::AbstractDouble) = exponent(x.hi)

baseof(::Type{Single{T}}) where T<:AbstractFloat = T
baseof(::Type{Double{T}}) where T<:AbstractFloat = T
baseof(::Type{Single{T}}) where T<:AbstractDouble = baseof(T)
baseof(::Type{Double{T}}) where T<:AbstractDouble = baseof(T)

realmax(::Type{T}) where T<:AbstractFloat = realmax(baseof(T))

maxintfloat(::Type{Single{T}}) where T<:AbstractFloat = Single(maxintfloat(T))
maxintfloat(::Type{Double{T}}) where T<:AbstractFloat = Double{T}(min(realmax(T),ldexp(1.0,precision(Double{T}))))

eps(::Type{Single{T}}) where T<:AbstractFloat = max(eps(zero(T)), eps(T)*eps(T))
eps(::Type{Double{T}}) where T<:AbstractFloat = max(eps(zero(T)), eps(T)*eps(T))

#typemin(::Type{TD}) where {TD <: AbstractDouble{T} where T <: AbstractFloat } = typemin(T)
#typemax(::Type{TD}) where {TD <: AbstractDouble{T} where T <: AbstractFloat } = typemax(T)

for D in [:Single, :Double]
  for op in [:typemin, :typemax]
    @eval $op(::Type{$D{T}}) where {T <: AbstractFloat} = $op(T)
  end
end

"""
eps(x::Double) returns the unit in last place, which may be very small since
there can be an arbitrary number of zero-bits between x.hi and x.lo.
"""
eps(x::Double) = eps(x.lo)
eps(x::Single) = eps(zero(x.hi))

"The minimum positive difference between possible rand() return values"
# NOTE rand(::Float64) in the current implementation always sets the last
# bit of the significand to zero. Hence eps() and not 0.5*eps()
rand_eps(::Type{Float64}) = eps(Float64)
rand_eps(::Type{Float32}) = eps(Float32)
rand_eps(::Type{Float16}) = eps(Float16)
rand_eps(::Type{Double{T}}) where T<:AbstractFloat = 2*rand_eps(T)*rand_eps(T)

function rand(rng::MersenneTwister, ::Type{Double{T}}) where {T<:AbstractFloat}
    normalize_double(rand(rng, T), rand_eps(T)*2*(rand(rng, T)-0.5))
end

function show(io::IO, x::AbstractDouble)
  if isnan(x)
    print(io, typeof(x), "(NaN)")
  elseif isinf(x)
    print(io, typeof(x), "(", x<0 ? "-" : "", "Inf)")
  else
    # correct number of decimal places, because convert to BigFloat sets precision.
    show(io, convert(BigFloat,x))
  end
end

# Not defined in this module, but should work anyway:
# zero, zeros, one, ones, sign...

end #module
