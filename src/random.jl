"The minimum positive difference between possible rand() return values"
# NOTE rand(::Float64) in the current implementation always sets the last
# bit of the significand to zero. Hence eps() and not 0.5*eps()
rand_eps(::Type{Float64}) = eps(Float64)
rand_eps(::Type{Float32}) = eps(Float32)
rand_eps(::Type{Float16}) = eps(Float16)
rand_eps{T}(::Type{Double{T}}) = 0.5*rand_eps(T)*rand_eps(T)

"""
A quick random number function, that does not fill all the bits of the
significand. (The spacing between possible random numbers is constant eps(2))
"""
function rand{T}(::Type{Double{T}})
    normalize_double(rand(T), rand_eps(T)*rand(T))
end

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

fullrand(::Type{Float32}) = Float32(fullrand(Float64));
fullrand(::Type{Float16}) = Float16(fullrand(Float64));

function fullrand{T}(::Type{Double{T}})
  hi = fullrand(T);
  lo = fullrand(T);
  return Double(hi,eps(hi)*(lo-0.5))
end

fullrand(::Type{Single}) = throw(InexactError())

function fullrand(::Type{BigFloat})
  x = fullrand(Float64)
  e = exponent(x)
  r = BigFloat(x)
  for i=1:Int(ceil((precision(BigFloat)-precision(Float64))/32))
    r = r + ldexp(Float64(rand(UInt32)),e-32*i)
  end
  return r
end
