using DoubleDoubles
using Base.Test

include("../src/fullrand.jl")

x = sqrt(2.0)
bx = big(x)
sx = Single(x)
dx = Double{Float64}(x)

y = 0.1
by = big(y)
sy = Single(y)
dy = Double{Float64}(y)

@test x == sx == dx
@test y == sy == dy

dxy = dx*dy
bxy = bx*by
@test sx*sy == dxy
@test x*y == Float64(dxy)
@test dxy == Double{Float64}(bxy)

@test x+y == Float64(dx+dy)
@test dx+dy == Double{Float64}(bx+by)

@test x-y == Float64(dx-dy)
@test dx-dy == Double{Float64}(bx-by)

@test x/y == Float64(dx/dy)
@test dx/dy == Double{Float64}(bx/by)

@test sqrt(y) == Float64(sqrt(dy))
@test sqrt(dy) == Double{Float64}(sqrt(by))

#@test rem(dxy,1.0) == Double(rem(bxy,1.0))


## New
@test Double{Float64}(pi) == Double{Float64}(3.141592653589793, 1.2246467991473532e-16)
@test Double{Float64}(3.5) == Double(3.5, 0.0)
@test Double{Float64}(3.5) == Double{Float64}(3.5, 0.0)

a = Double{Float64}(big"3.1")
@test a == Double{Float64}(3.1, -8.881784197001253e-17)

@test Single{Float64}(3) === Single(3.0)
@test Double{Float64}(3) == Double(3.0, 0.0)
@test Double{Float64}(big(3)) == Double(3.0, 0.0)

@test convert(Single{Float64}, 1) === Single(1.0)
@test convert(Single{Float32}, 1) === Single(1.0f0)
@test convert(Double{Float64}, 1) === Double(1.0, 0.0)
@test convert(Double{Float32}, 1) === Double(1.0f0, 0.0f0)

@test Single{Float32}(3) === Single{Float32}(3.0f0)
@test Double{Float32}(3) === Double{Float32}(3.0f0, 0.0f0)
@test Single{Float32}(BigFloat(3)) === Single{Float32}(3.0f0)
@test Double{Float32}(BigFloat(3)) === Double{Float32}(3.0f0, 0.0f0)

## Utility functions for test routine

testeps(x::AbstractFloat) = eps(x)
testeps(x::Integer) = one(x)

## Test routine for type interactions

function testtypes(A::Type, B::Type)
  a = rand(A)
  println("orig: ", a)
  if A<:Integer && B<:AbstractFloat && big(abs(a)) >= big(maxintfloat(B))
    a = mod(a,convert(A,maxintfloat(B)))
    println("1: ", a)
  end
  if A<:AbstractFloat && B<:Integer
    a = floor(min(maxintfloat(A),typemax(B))*a)
    println("2: ", a)
  end
  epsa = testeps(a)
  b = convert(B, a)
  epsb = testeps(b)
  a2 = convert(A, b)
  @test abs(big(a)-big(a2)) < 2*max(big(epsa), big(epsb))
end

## Test type interactions

intlist = [Int16, Int32, Int64, Int128, UInt16, UInt32, UInt64, UInt128, BigInt]
floatlist = [Float16, Float32, Float64, BigFloat]
singlelist = Array(DataType,0)
doublelist = Array(DataType,0)
for T = floatlist[1:end-1]
  for D = [Single{T}, Single{Single{T}}, Single{Single{Single{T}}},
          Single{Double{T}}, Single{Single{Double{T}}}, Single{Double{Double{T}}}]
    push!(singlelist, D)
  end
  for D = [Double{T}, Double{Double{T}}, Double{Double{Double{T}}}]
    push!(doublelist, D)
  end
end

typelist = [intlist; floatlist; singlelist; doublelist]

for A in typelist
  for B in typelist
    if (A <: AbstractDouble) || (B <: AbstractDouble)
      println(A, " - ", B)
      testtypes(A,B)
    end
  end
end
