"""
A slower random function that ensures that all bits of the significand are
randomly generated.

Any legal floating point number (excluding subnormals) in the range from 0 to
1-eps(T) might be returned (drawn from a uniform distribution).

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

function fullrand(::Type{Double{T}}) where T <: AbstractFloat
  hi = fullrand(T);
  lo = fullrand(T);
  return Double(hi,eps(hi)*(lo-1//2))
end

# Fullrand should use all bits, but for Single we make an exception.
fullrand(::Type{Single{T}}) where T<:AbstractFloat =
    Single(fullrand(T))

"fullrand(T,dims) returns an array of fullrand numbers of type T"
function fullrand(T::Type, dims::Dims)
  R = Array{T}(dims)
  for i in eachindex(R)
    R[i] = fullrand(T)
  end
  return R
end

fullrand(T::Type, d1::Integer, dims::Integer...) = fullrand(T, tuple(Int(d1), convert(Dims, dims)...))
