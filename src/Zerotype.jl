# We create this type just so that Single can have a .lo field.

"A type that stores no data, and holds the value zero."

immutable Zerotype <: AbstractFloat
end

promote_rule{T<:Number}(::Type{Zerotype}, ::Type{T}) = T

convert(::Type{Zerotype}, ::Zerotype) = Zerotype()
convert{T<:Number}(::Type{T}, ::Zerotype) = zero(T)
convert{T<:Number}(::Type{Zerotype}, x::T) = x==zero(T) ? Zerotype() : throw(InexactError)

# (Some of these rules do not handle Inf and NaN according to the IEEE spec.)

+(x::Number,::Zerotype) = x
+(::Zerotype,x::Number) = x
+(x::Number,::Zerotype) = x
+(::Zerotype,::Zerotype) = Zerotype()

-(x::Number,::Zerotype) = x
-(::Zerotype,x::Number) = -x
-(::Zerotype,::Zerotype) = Zerotype()
-(::Zerotype) = Zerotype()

*(::Number,::Zerotype) = Zerotype()
*(::Zerotype,::Number) = Zerotype()
*(::Zerotype,::Zerotype) = Zerotype()

/(::Zerotype,::Number) = Zerotype()

<(::Zerotype,::Zerotype) = false
>(::Zerotype,::Zerotype) = false
>=(::Zerotype,::Zerotype) = true
<=(::Zerotype,::Zerotype) = true

show(io::IO, ::Zerotype) = print(io, "0̸")
