"Combination of maxintfloat(AbstractFloat) and typemax(Integer)"
tmax(::Type{T}) where {T <: AbstractFloat} = maxintfloat(T)
tmax(::Type{T}) where {T <: Integer} = typemax(T)

"a version of `abs` that will never return negative numbers"
aabs(x::Integer) = x==typemin(x) ? zero(x) : abs(x)

"""
A random number useful for testing operations on a type.
This number does not follow any particular distribution, but should
cover a lot of possible bit-patterns, while excluding `Inf`, `NaN`, and
subnormals.
"""
testrand(T) = rand(T)
function testrand(::Type{T}) where T <: Base.Math.IEEEFloat
  u = rand(Base.fpinttype(T))
  if u & Base.exponent_mask(T) == Base.exponent_mask(T)
    u = u & (Base.exponent_one(T) | ~Base.exponent_mask(T))
  end
  return reinterpret(T, u)
end
testrand(::Type{Single{T}}) where {T <: AbstractFloat} = Single(testrand(T))
function testrand(::Type{Double{T}}) where {T <: AbstractFloat}
  hi = testrand(T);
  lo = eps(hi)*fullrand(T);
  return Double(hi,lo)
end

function testrand(::Type{BigInt})
  r = BigInt(0)
  for i=1:8
    r = r<<64+rand(Int64)
  end
  return r
end

testrand(::Type{BigFloat}) = ldexp(BigFloat(testrand(BigInt)), rand(Int64)%4611686018427387904)

"Minimum of the max int supported by two types"
minmaxint(::Type{T1},::Type{T2}) where {T2 <: AbstractFloat} where {T1 <: Integer} =
    min(big(typemax(T1)),big(maxintfloat(T2)))
minmaxint(::Type{BigInt},::Type{T2}) where {T2 <: AbstractFloat} =
    big(maxintfloat(T2))

"Reflect that the true precision of a Single is less than `precision` indicates"
trueprecision(T::Type) = precision(T)
trueprecision(::Type{Single{T}}) where {T <: AbstractFloat} = trueprecision(T)

"""
Given two types, `testrand` returns a random number of the first type,
such that it can be converted into the second type without loss of accuracy.
(Useful for testing conversion.)
"""
testrand(::Type{T}, ::Type{T}) where {T<:Any} = testrand(T)

testrand(::Type{T1}, ::Type{T2}) where {T2 <: Unsigned} where {T1 <: Unsigned} = testrand(T2) % T1
testrand(::Type{T1}, ::Type{T2}) where {T2 <: Integer} where {T1 <: Unsigned} = testrand(T2) % T1
testrand(::Type{T1}, ::Type{T2}) where {T2 <: Unsigned} where {T1 <: Integer} = aabs(testrand(T2) % T1)
testrand(::Type{T1}, ::Type{T2}) where {T2 <: Integer} where {T1 <: Integer} = testrand(T2) % T1

function testrand(::Type{T1}, ::Type{T2}) where {T2 <: AbstractFloat} where {T1 <: Integer}
  return testrand(T1)%convert(T1,minmaxint(T1,T2))
end

function testrand(::Type{T1}, ::Type{T2}) where {T2 <: Integer} where {T1 <: AbstractFloat}
  return convert(T1, testrand(T2, T1))
end

testrand(::Type{BigFloat}, ::Type{T}) where {T <: AbstractFloat} = fullrand(T)
testrand(::Type{T}, ::Type{BigFloat}) where {T <: AbstractFloat} = testrand(T)
testrand(::Type{T}, ::Type{BigInt}) where {T <: AbstractFloat} =
    testrand(T, Int64)

function testrand(::Type{T1}, ::Type{T2}) where {T2 <: AbstractFloat} where {T1 <: AbstractFloat}
  if trueprecision(T1) < trueprecision(T2)
    return fullrand(T1)
  else
    return convert(T1, fullrand(T2))
  end
end

"testrand(T1,dims) returns an array of fullrand numbers of type T"
function testrand(T1::Type, dims::Dims)
  R = Array{T}(dims)
  for i in eachindex(R)
    R[i] = testrand(T1)
  end
  return R
end

testrand(T1::Type, T2::Type, d1::Integer, dims::Integer...) = testrand(T, tuple(Int(d1), convert(Dims, dims)...))

"testrand(T1,T2,dims) returns an array of fullrand numbers of type T"
function testrand(T1::Type, T2::Type, dims::Dims)
  R = Array{T}(dims)
  for i in eachindex(R)
    R[i] = testrand(T1, T2)
  end
  return R
end

testrand(T1::Type, T2::Type, d1::Integer, dims::Integer...) = testrand(T, tuple(Int(d1), convert(Dims, dims)...))
