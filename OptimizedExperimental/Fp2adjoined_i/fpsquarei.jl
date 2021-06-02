import Primes: isprime
import Base: +,-,*,/,^,==,inv,show

# Field has only one attribute, p
# Implicitly Fp^2 initialized as Fp/<x^2+1>
# Requires p = 3 (mod 4) for polynomial to be irredusible
struct Fp2{T <: Integer}
    p::T
    #For calculating cuberoots. k = (p^2+2)/9
    k::T
end

# FieldElements point to their parent field
# Consists of integers a,b representing the integer a + b*i
struct Fp2Elem{T <: Integer}
    a::T
    b::T
    field::Fp2{T}
end

# Field-Constructor.
function Fp2(p::T) where T <: Integer
    !isprime(p) && error("p must prime")
    mod(p, 4) != 3 && error("p must be 3 (mod 4)")
    if mod(p, 9) == 4 || mod(p, 9) == 5
        return Fp2{T}(p, ((p*p)+2)÷9)
    else
        error("p must be +-5 (mod 9)")
    end
end

# Element-Constructor
# ASSUMES a,b are both between 0 and p for speed.
# Arithmetic functions will make sure of this
function (F::Fp2{T})(a::T,b::T) where T <: Integer
    return Fp2Elem{T}(a, b, F)
end

############################
#
# Printing
#
############################

function show(io::IO, F::Fp2{T}) where {T <: Integer}
    print(io, "F_$(F.p)[X]/<X^2+1>")
end

function show(io::IO, F::Fp2Elem{T}) where {T <: Integer}
    print(io, "$(F.a) + $(F.b)*i")
end

############################
#
# Some helper functions
#
############################

function zero(F::Fp2{T}) where T <: Integer
    return F(T(0),T(0))
end

function one(F::Fp2{T}) where T <: Integer
    return F(T(1),T(0))
end

function iszero(x::Fp2Elem{T}) where T <: Integer
    return x.a == 0 && x.b == 0
end

function isone(x::Fp2Elem{T}) where T <: Integer
    return x.a == 1 && x.b == 0
end

function ==(x::Fp2Elem{T}, y::Fp2Elem{T}) where T <: Integer
    return (x.a == y.a) && (x.b == y.b)
end

#Assumes a is a cube
function cuberoot(x::Fp2Elem{T}) where T <: Integer
    return x^x.field.k
end

############################
#
# Basic operators
#
############################
# All operators trusts that elements are from the same ring
# Additive inverse
function -(x::Fp2Elem{T}) where T <: Integer
    a_new = mod(-x.a, x.field.p)
    b_new = mod(-x.b, x.field.p)
    return x.field(a_new, b_new)
end

# Addition
function +(x::Fp2Elem{T}, y::Fp2Elem{T}) where T <: Integer
    a_new = mod((x.a + y.a), x.field.p)
    b_new = mod((x.b + y.b), x.field.p)
    return x.field(a_new, b_new)
end

# Subtraction
function -(x::Fp2Elem{T}, y::Fp2Elem{T}) where T <: Integer
    return x+(-y)
end

# Multiplication
function *(x::Fp2Elem{T}, y::Fp2Elem{T}) where T <: Integer
    a_new = mod((x.a*y.a - x.b*y.b), x.field.p)
    b_new = mod((x.a*y.b + y.a*x.b), x.field.p)
    return x.field(a_new, b_new)
end

# Scalar multiplication on left
function *(x::Integer, y::Fp2Elem{T}) where T <: Integer
    return y.field(mod(T(x), y.field.p),T(0))*y
end

# Scalar multiplication on right
function *(x::Fp2Elem{T}, y::Integer) where T <: Integer
    return x*x.field(mod(T(y), x.field.p),T(0))
end

# Multiplicative inverse
function inv(x::Fp2Elem{T}) where T <: Integer
    x.a == 0 && x.b == 0 && error("0 is not invertible...")
    #x^-1 = x.a*(x.a^2+x.b^2)^-1 + (-x.b)(x.a^2+x.b^2)^-1 * i
    d = invmod((x.a*x.a + x.b*x.b),x.field.p)
    a_new = mod((x.a*d), x.field.p)
    b_new = mod((-x.b*d), x.field.p)
    return x.field(a_new, b_new)
end

# Exponensiation
# Uses the simple square and multiply algorithm
function ^(x::Fp2Elem{T}, y::Integer) where T <: Integer
    if iszero(x) || y == 1
        return x
    end
    if y == 0
        return one(x.field)
    end
    if y < 0
        x = inv(x)
        y = -y
    end
    bit = T(1) << (ndigits(y, base = 2) - 1)
    z = x
    bit >>= 1
    while bit != 0
        z = z*z
        if (bit & y) != 0
            z *= x
        end
        bit >>=1
    end
    return z
end
