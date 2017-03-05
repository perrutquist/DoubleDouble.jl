module DoubleDoubles

export Double, Single, DoubleDouble, QuadDouble, OctaDouble
import Base:
    convert,
    *, +, -, /, sqrt, <, >, <=, >=, ldexp,
    rem, abs, rand, promote_rule,
    show, big,
    precision

include("Zerotype.jl")

abstract AbstractDouble{T} <: AbstractFloat

"A Single is a Double where the lower part is exactly zero."
immutable Single{T<:AbstractFloat} <: AbstractDouble{T}
    hi::T
    lo::Zerotype
end
Single(hi::AbstractFloat) = Single(hi,Zerotype())

"Double{T} uses two instances of type T to represent a floating point number."
# In a Double, hi uses the full mantissa, and abs(lo) <= 0.5eps(hi)
# This must be enforced by the method that creates the double.
immutable Double{T<:AbstractFloat} <: AbstractDouble{T}
    hi::T
    lo::T
end

function Double(x::Real)
  warn("Double(x) is depreciated. Use Double{Float64}(x) or DoubleDouble(x) to convert x to Double{Float64}.")
  convert(Double{Float64}, x)
end

# A Double with a Zerotype as lo is by definition a Single.
Double(hi::AbstractFloat, ::Zerotype) = Single(hi)

# "Normalise" Doubles to ensure abs(lo) <= 0.5eps(hi)
# assumes abs(u) > abs(v): if not, use Single + Single
function normalize_double{T}(u::T, v::T)
    w = u + v
    Double(w, (u-w) + v)
end

# This is a lower bound on the precision. If the binary representation of a
# number has a lot of zeros in the middle, then higher precision is achieved.
# For example, 1 + 2^-300 can be represented exactly by a Double{Float64}
precision{T}(::Type{Double{T}}) = 2*precision(T)

# The precision of Single is defined the same as the precision of the
# correspoinding Double. The digits "stored" in the always-zero .lo
# field count! (This is why Single contructors may throw InexactError.)
precision{T}(::Type{Single{T}}) = 2*precision(T)

# ldexp can be done without renormalization
ldexp(x::Double, n::Integer) = Double(ldexp(x.hi,n), ldexp(x.lo,n))

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
@generated halfprec_constant{T}(::Type{T}) = convert(T,exp2(ceil(precision(T)/2)))

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

## conversion

convert{T<:AbstractFloat}(::Type{Single{T}}, x::T) = Single(x)
convert{T<:AbstractFloat}(::Type{Double{T}}, x::T) = Double(x, zero(T))

convert{T<:AbstractFloat}(::Type{Double{T}}, x::Single{T}) = Double(x.hi, zero(T))
convert{T<:AbstractFloat}(::Type{Single{T}}, x::Double{T}) = x.lo == zero(T) ? Single(x.hi) : throw(InexactError())

convert{T<:AbstractFloat}(::Type{T}, x::AbstractDouble{T}) = x.hi

convert{T<:AbstractFloat}(::Type{Single{T}}, x::Single{T}) = x # needed because Single <: AbstractFloat
convert{T<:AbstractFloat}(::Type{Double{T}}, x::Double{T}) = x # needed because Double <: AbstractFloat

function convert{T<:AbstractFloat}(::Type{Single{T}}, x::Real)
  z = convert(T,x)
  z == x || throw(InexactError())
  Single(z)
end

function convert{T<:AbstractFloat}(::Type{Single{T}}, x::Rational)
  z = convert(T,x)
  z == x || throw(InexactError())
  Single(z)
end

convert{T1<:AbstractFloat, T2<:AbstractFloat}(::Type{T1}, x::AbstractDouble{T2}) =
  convert(T1, x.hi) + (precision(T1) >= precision(T2) ? convert(T1, x.lo) : 0)

function convert{T<:AbstractFloat}(::Type{Double{T}}, x::Real)
  z = convert(T, x)
  Double(z, convert(T, x-z))
end

function convert{T<:AbstractFloat}(::Type{BigFloat}, x::AbstractDouble{T})
  setprecision(BigFloat, 2*precision(T)) do
    convert(BigFloat, x.hi) + convert(BigFloat, x.lo)
  end
end

convert{T<:AbstractFloat}(::Type{Double{T}}, x::Rational) = convert(Double{T}, x.num)/convert(Double{T}, x.den)

## promotion
promote_rule{T<:AbstractFloat}(::Type{Double{T}}, ::Type{Int64}) = Double{T}
promote_rule{T<:AbstractFloat}(::Type{Single{T}}, ::Type{Int64}) = Double{T}

promote_rule{T<:AbstractFloat}(::Type{Single{T}}, ::Type{T}) = Single{T}
promote_rule{T<:AbstractFloat}(::Type{Double{T}}, ::Type{T}) = Double{T}
promote_rule{T<:AbstractFloat}(::Type{Double{T}}, ::Type{Single{T}}) = Double{T}

promote_rule{s,T<:AbstractFloat}(::Type{Irrational{s}}, ::Type{Single{T}}) = Double{T}

# We'll assume that AbstractDouble is used for speed and BigFloat for precision
# The construct Double{BigFloat} should be avoided in favour of selecting
# higher bigfloat precision.
promote_rule{T<:AbstractFloat}(::Type{AbstractDouble{T}}, ::Type{BigFloat}) = BigFloat

# Double(::Type{T}, x::Real) converts x to Double{T}, computing both the hi and lo
# parts. There's no constructor that explicitly puts a zero in the lo.
# Use Double(x, zero(x)) if Single(x) cannot be used.

Double{T}(::Type{T}, x::Real) = convert(Double{T}, x)

# Macro will do the BigFloat conversion at compile-time.
@generated function doublesym{T<:AbstractFloat,sym}(::Type{Double{T}}, ::Irrational{sym})
  Double(T, setprecision(()->BigFloat(eval(sym)),BigFloat,2*precision(T)+32))
end

convert{T<:AbstractFloat}(::Type{Double{T}}, sym::Irrational) = doublesym(Double{T}, sym)

# Type aliases

typealias DoubleDouble Double{Float64}
typealias Quad{T} Double{Double{T}}
typealias QuadDouble Quad{Float64}
typealias Octa{T} Double{Quad{T}}
typealias OctaDouble Octa{Float64}

# Comparisons

for op in [:<,:>]
  @eval $op{T}(x::AbstractDouble{T},y::AbstractDouble{T})=$op(x.hi,y.hi)?true:$op(x.lo,y.lo)
end
for op in [:<=,:>=]
  @eval $op{T}(x::AbstractDouble{T},y::AbstractDouble{T})=$op(x.hi,y.hi)?$op(x.lo,y.lo):false
end

# TODO eliminate branches
# add12
function +{T}(x::Single{T},y::Single{T})
    abs(x.hi) > abs(y.hi) ? normalize_double(x.hi, y.hi) : normalize_double(y.hi, x.hi)
end

# Dekker add2
function +{T}(x::AbstractDouble{T}, y::AbstractDouble{T})
    r = x.hi + y.hi
    s = abs(x.hi) > abs(y.hi) ? (((x.hi - r) + y.hi) + y.lo) + x.lo : (((y.hi - r) + x.hi) + x.lo) + y.lo
    normalize_double(r, s)
end

-{T<:AbstractFloat}(x::AbstractDouble{T}) = Double(-x.hi, -x.lo)

function -{T}(x::Single{T},y::Single{T})
    abs(x.hi) > abs(y.hi) ? normalize_double(x.hi, -y.hi) : normalize_double(-y.hi, x.hi)
end

function -{T}(x::AbstractDouble{T}, y::AbstractDouble{T})
    r = x.hi - y.hi
    s = abs(x.hi) > abs(y.hi) ? (((x.hi - r) - y.hi) - y.lo) + x.lo : (((-y.hi - r) + x.hi) + x.lo) - y.lo
    normalize_double(r, s)
end

# TODO FMA version
# Dekker mul12
function *{T}(x::Single{T},y::Single{T})
    hx,lx = splitprec(x.hi)
    hy,ly = splitprec(y.hi)
    z = x.hi*y.hi
    normalize_double(z, ((hx*hy-z) + hx*ly + lx*hy) + lx*ly)
end

# Dekker mul2
function *{T}(x::AbstractDouble{T}, y::AbstractDouble{T})
    c = Single(x.hi) * Single(y.hi)
    cc = (x.hi * y.lo + x.lo * y.hi) + c.lo
    normalize_double(c.hi, cc)
end

# Dekker div2
function /{T}(x::AbstractDouble{T}, y::AbstractDouble{T})
    c = x.hi / y.hi
    u = Single(c) * Single(y.hi)
    cc = ((((x.hi - u.hi) - u.lo) + x.lo) - c*y.lo)/y.hi
    normalize_double(c, cc)
end

# Dekker sqrt2
function sqrt{T}(x::AbstractDouble{T})
    if x.hi <= 0
        throw(DomainError("sqrt will only return a complex result if called with a complex argument."))
    end
    c = sqrt(x.hi)
    u = Single(c)*Single(c)
    cc = (((x.hi - u.hi) - u.lo) + x.lo)*convert(T,0.5)/c
    Double(c, cc)
end

rem{T}(x::Double{T},d::Real) = normalize_double(rem(x.hi,d), rem(x.lo,d))
abs{T}(x::Double{T}) = x.hi>0?x:-x

include("random.jl")

function show{T}(io::IO, x::AbstractDouble{T})
  print(io, convert(BigFloat,x))
end

end #module
