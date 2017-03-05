"""
A slower random function that ensures that all bits of the significand are
randomly generated.

Any legal floating point number (excluding subnormals) in the range from 0 to
1-eps() might be returned (drawn from a uniform distribution).

The spacing beteeen possible random numbers near x is eps(x).

If rand() is used in testing arithmetic operations, then trailing zeros
in the significand might mask potential problems. With fullrand, this
problem is avoided.
"""

function fullrand(::Type{Float64})
  rnd = rand(UInt64)
  lz = leading_zeros(rnd)
  while rnd==0
    rnd = rand(UInt64)
    lz += leading_zeros(rnd)
    if lz>=1022
      return 0.0 #If we get here then rand(Uint64) is probably broken...
    end
  end
  return reinterpret(Float64, ((1022-lz)<<52) | (rand(UInt64)&(2^52-1)))
end

fullrand(::Type{Float32}) = Float32(fullrand(Float64))
fullrand(::Type{Float16}) = Float16(fullrand(Float64))

function fullrand{T}(::Type{Double{T}})
  hi = fullrand(T);
  lo = fullrand(T);
  return Double(hi,eps(hi)*(lo-0.5))
end

# Really, fullrand should use all bits, so we should have:
# fullrand(::Type{Single}) = throw(InexactError())
# But since we're only using it for testing, it makes sense to have
# a fullrand() that returns a Single
fullrand{T}(::Type{Single{T}}) = Single(fullrand(T))

function fullrand(::Type{BigFloat})
  x = fullrand(Float64)
  e = exponent(x)
  r = BigFloat(x)
  for i=1:Int(ceil((precision(BigFloat)-precision(Float64))/32))
    r = r + ldexp(Float64(rand(UInt32)),e-32*i)
  end
  return r
end

"fullrand() for integer types is an alias for rand()"
fullrand{T<:Integer}(::Type{T}) = rand(T)

# fullrand for BigInt is a misnomer. BigInt has infinite precision,
# but obviously we cannot generate an infinite number of random bits.
# The goal here is to generate something that is useful for testing.
"""
fullrand(BigInt) generates some very large integers, but that will still not
cause overflow in a Float64.
"""
function fullrand(::Type{BigInt})
  r = BigInt(0)
  for i=1:8
    r = r<<64+rand(Int64)
  end
  return r
end

"fullrand(T,dims) returns an array of fullrand numbers of type T"
function fullrand(T::Type, dims::Dims)
  R = Array{T}(dims)
  for i in eachindex(R)
    R[i] = fullrand(T)
  end
  return R
end

fullrand(T::Type, d1::Integer, dims::Integer...) = fullrand(T, tuple(Int(d1), convert(Dims, dims)...))
