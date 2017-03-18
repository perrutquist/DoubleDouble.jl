# module Zeros

# rename Zerotype -> Zero and make this a module!

import Base: +, -, *, /, <, >, <=, >=,
     ldexp, copysign, flipsign, sign, round, floor, ceil, trunc,
     promote_rule, convert, show, significand, isodd, iseven

#export Zero ifzero

# Preformance comparison: Naïve implementation of axpy!()
# axpy!(a,X,Y) = Y .= a.*X .+ Y
# r = rand(Float64, 1000_000);
# @time axpy(0.0, r, r)
# @time axpy(Zerotype(), r, r)

 "A type that stores no data, and holds the value zero."
immutable Zerotype <: Real
end

promote_rule{T<:Number}(::Type{Zerotype}, ::Type{T}) = T

convert(::Type{Zerotype}, ::Zerotype) = Zerotype()
convert{T<:Real}(::Type{T}, ::Zerotype) = zero(T)
convert{T<:Real}(::Type{Zerotype}, x::T) = x==zero(T) ? Zerotype() : throw(InexactError)

# (Some of these rules do not handle Inf and NaN according to the IEEE spec.)

+(x::Number,::Zerotype) = x
+(::Zerotype,x::Number) = x
+(::Zerotype,::Zerotype) = Zerotype()

-(x::Number,::Zerotype) = x
-(::Zerotype,x::Number) = -x
-(::Zerotype,::Zerotype) = Zerotype()

*(::Number,::Zerotype) = Zerotype()
*(::Zerotype,::Number) = Zerotype()
*(::Zerotype,::Zerotype) = Zerotype()

/(::Zerotype,::Number) = Zerotype()

<(::Zerotype,::Zerotype) = false
>(::Zerotype,::Zerotype) = false
>=(::Zerotype,::Zerotype) = true
<=(::Zerotype,::Zerotype) = true
# == and != work by default methods

ldexp(::Zerotype, ::Integer) = Zerotype()
copysign(::Zerotype,::Real) = Zerotype()
flipsign(::Zerotype,::Real) = Zerotype()

# Alerady working due to default definitions:
# abs, isinf, isnan, isinteger, isreal, isimag,
# signed, unsigned, float, big, complex

for op in [:(-), :sign, :round, :floor, :ceil, :trunc, :significand]
  @eval $op(::Zerotype) = Zerotype()
end

isodd(::Zerotype) = false
iseven(::Zerotype) = true

show(io::IO, ::Zerotype) = print(io, "0̸")

"Convert to Zerotype() if zero. (Use immediately before calling a function.)"
# This function is intentionally not type-stable.
testzero(x::Real) = x==0 ? Zerotype() : x

#end #module
