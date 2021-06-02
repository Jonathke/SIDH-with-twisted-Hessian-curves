module FieldArithmetic

export Ring, RingElem, CRing, CRingElem, Field, FieldElem, Zmod, ZmodElem, PolyRing, PolyRingElem,
QuotientRing, QuotientRingElem, Fpsquare, order, squareroot, cuberoot, cuberootOfOne, allcuberoots,
allElements

##########################################
#         For overwriting
##########################################

import Base: +,-,*,/,^,==,inv,show,zeros,gcd,gcdx,iszero,isone,zero,one
import Primes: isprime

##########################################
# Abstract types for defining domains
##########################################

abstract type Ring end

abstract type RingElem end

abstract type CRing <: Ring end

abstract type CRingElem <: RingElem end

##########################################
#        Generic helper functions
##########################################

function sameRing(x::RingElem, y::RingElem)
    return x.ring == y.ring
end

function zeros(R::T, len::Integer) where T <: Ring
    return [zero(R) for i in 1:len]
end

##########################################
#
#         Ring of Integers mod n
#
##########################################

struct Zmod{T <: Integer} <: CRing
    n::T
end

struct ZmodElem{T <: Integer} <: CRingElem
    d::T
    ring::Zmod{T}
end

# Ring-Constructor
function Zmod(n::T) where T <: Integer
    return Zmod{T}(n)
end

# Element-Constructor
function (R::Zmod{T})(a::Integer) where T <: Integer
    a = T(a)
    n = R.n
    d = a%n
    if d < 0
        d += n
    end
    return ZmodElem{T}(d, R)
end

###### Printing ######
function show(io::IO, R::Zmod{T}) where T <: Integer
    print(io, "Z_$(R.n)")
end

function show(io::IO, x::ZmodElem{T}) where T <: Integer
    print(io, "$(x.d)")
end

###### Helper functions ######

function zero(R::Zmod{T}) where T <: Integer
    return R(T(0))
end

function one(R::Zmod{T}) where T <: Integer
    return R(T(1))
end

function iszero(x::ZmodElem{T}) where T <: Integer
    return x.d == 0
end

function isone(x::ZmodElem{T}) where T <: Integer
    return x.d == 1
end

function order(R::Zmod{T}) where T <: Integer
    return R.n
end

function allElements(R::Zmod{T}) where T <: Integer
    return [R(i) for i in 0:R.n-1]
end

###### Basic operators #######
function -(x::ZmodElem{T}) where T <: Integer
    if x.d == 0
        return x
    else
        R = x.ring
        return R(R.n - x.d)
    end
end

function +(x::ZmodElem{T}, y::ZmodElem{T}) where T <: Integer
    if !sameRing(x, y)
        error("Cannot add elements of different rings")
    end
    z = x.d + y.d
    R = x.ring
    return R(z)
end

function -(x::ZmodElem{T}, y::ZmodElem{T}) where T <: Integer
    return x+(-y)
end

function *(x::ZmodElem{T}, y::ZmodElem{T}) where T <: Integer
    if !sameRing(x, y)
        error("Cannot multiply elements of different rings")
    end
    z = x.d * y.d
    R = x.ring
    return R(z)
end

function *(x::Integer, y::ZmodElem{T}) where T <: Integer
    return y.ring(x)*y
end

function *(x::ZmodElem{T}, y::Integer) where T <: Integer
    return x.ring(y)*x
end

function inv(x::ZmodElem{T}) where T <: Integer
    R = x.ring
    g, a, b = gcdx(x.d, R.n)
    g != 1 && error("Element $x is not invertible in $(x.ring)")
    return R(a)
end

#exponensiation uses the simple square and multiply algorithm
function ^(x::ZmodElem{T}, y::Integer) where T <: Integer
    if x.d == 0 || y == 1
        return x
    end
    R = x.ring
    if y == 0
        return R(T(1))
    end
    if y < 0
        x = inv(x)
        y = -y
    end
    # Should include some reduction of y modulo order if order is known
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

function ==(x::ZmodElem{T}, y::ZmodElem{T}) where T <: Integer
    return sameRing(x, y) && (x.d == y.d)
end

function ==(R1::Zmod{T}, R2::Zmod{T}) where T <: Integer
    return R1.n == R2.n
end



##########################################
#
#         Polynomial Ring
#
##########################################

struct PolyRing{T <: CRing} <: CRing
    parent::T
    variable::String
end

struct PolyRingElem{T <: CRingElem, S <: CRing} <: CRingElem
    d::Array{T}
    ring::PolyRing{S}
end

# Ring-Constructor
function PolyRing(R::T,var::String) where T <: CRing
    return PolyRing{T}(R,var)
end

function PolyRing(R::T) where T <: CRing
    return PolyRing{T}(R,"X")
end

# Element-Constructor
function (R::PolyRing{S})(a::Array{T}) where {T <: CRingElem, S <: CRing}
    while length(a) > 0
        if !(iszero(a[length(a)]))
            return PolyRingElem{T,S}(a, R)
        end
        pop!(a)
    end
    return PolyRingElem{T,S}([zero(R.parent)], R)
end

function (R::PolyRing{Zmod{T}})(a::Array{S}) where {T <: Integer, S <: Integer}
    Z = R.parent
    a = [Z(i) for i in a]
    return R(a)
end

###### Printing ######
function show(io::IO, R::PolyRing{T}) where T <: CRing
    print(io, "$(R.parent)[$(R.variable)]")
end

function show(io::IO, x::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    s = "$(x.d[1])"
    l = length(x.d)
    R = x.ring
    if l > 1
        s *= " + $(x.d[2])*$(R.variable)"
        if l > 2
            for i in 3:l
                s *= " + $(x.d[i])*$(R.variable)^$(i-1)"
            end
        end
    end
    print(io, s)
end

###### Useful helper functions #######
function zero(R::PolyRing{T}) where T <: CRing
    return R([zero(R.parent)])
end

function one(R::PolyRing{T}) where T <: CRing
    return R([one(R.parent)])
end

function iszero(a::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    return a==zero(a.ring)
end

function isone(a::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    return a==one(a.ring)
end

function degree(a::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    if iszero(a)
        return -1
    else
        return length(a.d)-1
    end
end

function LC(a::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    return a.d[length(a.d)]
end

function monic(a::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    x = inv(LC(a))
    return a*x
end

###### Basic operators #######
function -(a::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    d = deepcopy(a.d)
    l = length(d)
    R = a.ring
    z = zeros(R.parent, l)
    for i in 1:l
        z[i] = -d[i]
    end
    return R(z)
end

function +(a::PolyRingElem{T,S}, b::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    if !sameRing(a, b)
        error("Cannot add elements from different rings")
    end
    R = a.ring
    d1 = deepcopy(a.d)
    d2 = deepcopy(b.d)
    lena = length(d1)
    lenb = length(d2)
    lenz = max(lena,lenb)
    lenmin = min(lena,lenb)
    z = zeros(R.parent, lenz)
    for i in 1:lenmin
        z[i] = d1[i] + d2[i]
    end
    if lena > lenmin
        for i in lenmin+1:lena
            z[i] = d1[i]
        end
    end
    if lenb > lenmin
        for i in lenmin+1:lenb
            z[i] = d2[i]
        end
    end
    return R(z)
end

function -(a::PolyRingElem{T,S}, b::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    return a + (-b)
end

#Might be subject to more effective implementations later...
function *(a::PolyRingElem{T,S}, b::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    if !sameRing(a, b)
        error("Cannot multiply elements from different rings")
    end
    R = a.ring
    d1 = deepcopy(a.d)
    d2 = deepcopy(b.d)
    lena = length(d1)
    lenb = length(d2)
    if lena == 0 || lenb == 0
        return R([zero(R.parent)])
    end
    lenz = lena + lenb - 1
    z = zeros(R.parent, lenz)
    k = lenz-lena
    for i in 1:k
        push!(d1,zero(R.parent))
    end
    k = lenz-lenb
    for i in 1:k
        push!(d2,zero(R.parent))
    end
    for i in 1:lenz
        s = zero(R.parent)
        for j in 1:i
            s += d1[j] * d2[i+1-j]
        end
        z[i] = s
    end
    return R(z)
end

function *(a::PolyRingElem{T,S}, b::T) where {T <: CRingElem, S <: CRing}
    R = a.ring
    z = deepcopy(a.d)
    z = [k*b for k in z]
    return R(z)
end

function *(b::PolyRingElem{T,S}, a::T) where {T <: CRingElem, S <: CRing}
    R = b.ring
    z = deepcopy(b.d)
    z = [k*a for k in z]
    return R(z)
end

function *(n::U, a::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing, U <:Integer}
    if n < 0
        a = -a
        n = -n
    end
    if n == 1
        return a
    end
    Rpoly = a.ring
    R = Rpoly.parent
    if n == 0
        return zero(Rpoly)
    end
    bit = U(1) << (ndigits(n, base = 2) - 1)
    z = a
    bit >>= 1
    while bit != 0
        z = z+z
        if (bit & n) != 0
            z += a
        end
        bit >>=1
    end
    return z
end

function *(a::PolyRingElem{T,S}, n::U) where {T <: CRingElem, S <: CRing, U <:Integer}
    return n*a
end

#Again, just using square and multiply algorithm
function ^(a::PolyRingElem{T,S}, y::U) where {T <: CRingElem, S <: CRing, U <: Integer}
    if y < 0
        #Will implement invert function latert
        error("Not yet supported")
    end
    if y == 1
        return x
    end
    Rpoly = a.ring
    R = Rpoly.parent
    if y == 0
        return one(Rpoly)
    end
    # Should include some reduction of y modulo order if order is known
    bit = U(1) << (ndigits(y, base = 2) - 1)
    z = a
    bit >>= 1
    while bit != 0
        z = z*z
        if (bit & y) != 0
            z *= a
        end
        bit >>=1
    end
    return z
end

function ==(a::PolyRingElem{T,S}, b::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    return sameRing(a, b) && (a.d == b.d)
end

####### GCD (Very much needed for constructing quotient rings) #######
function gcd(a::PolyRingElem{T,S}, b::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    if !sameRing(a, b)
        error("Elements must be elements of same ring")
    end
    if degree(a) < degree(b)
        return gcd(b,a)
    end
    if iszero(a) && iszero(b)
        return zero(a.ring)
    end
    while degree(b) > -1
        q,r = divmod(a,b)
        a = b
        b = r
    end
    return monic(a)
end

function gcdx(a::PolyRingElem{T,S}, b::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    if !sameRing(a, b)
        error("Elements must be elements of same ring")
    end
    if degree(a) < degree(b)
        (a, y1, x1) = gcdx(b,a)
        return (a, x1, y1)
    end
    Rpol = a.ring
    if iszero(b)
        factor = inv(LC(a))
        return (a*factor, one(Rpol)*factor, zero(Rpol))
    end
    x0, x1 = zero(Rpol), one(Rpol)
    y1, y0 = zero(Rpol), one(Rpol)
    while degree(b) > -1
        q, r = divmod(a,b)
        x = x1 - q*x0
        y = y1 - q*y0
        a, b = b, r
        x0, x1, y0, y1 = x, x0, y, y0
    end
    factor = inv(LC(a))
    return (a*factor, x1*factor, y1*factor)
end

#gcdx subroutine
function divmod(a::PolyRingElem{T,S}, b::PolyRingElem{T,S}) where {T <: CRingElem, S <: CRing}
    if !sameRing(a, b)
        error("Elements must be elements of same ring")
    end
    Rpol = a.ring
    R = Rpol.parent
    lcdiv = LC(b)
    quotient = zero(Rpol)
    remainder = a
    while degree(remainder) >= degree(b)
        z = zeros(R, degree(remainder)-degree(b))
        push!(z,LC(remainder)*(inv(lcdiv)))
        z = Rpol(z)
        quotient += z
        remainder -= z*b
    end
    return quotient, remainder
end


##########################################
#
#         Quotient Ring
#
##########################################


struct QuotientRing{T <: CRing, S <: CRingElem} <: CRing
    parent::T
    modulo::S
end

struct QuotientRingElem{T <: CRing, S <: CRingElem} <: CRingElem
    d::S
    ring::QuotientRing{T,S}
end

# Ring-Constructor
function QuotientRing(R::T, m::S) where {T <: CRing, S <: CRingElem}
    return QuotientRing{T,S}(R,m)
end

# Element-Constructor
function (R::QuotientRing{T,S})(a::S) where {T <: CRing, S <: CRingElem}
    m = R.modulo
    q,r = divmod(a,m)
    return QuotientRingElem{T,S}(r, R)
end

function (R::QuotientRing{PolyRing{T},PolyRingElem{S,T}})(a::Array{S}) where {T <: CRing, S <: CRingElem}
    a = R.parent(a)
    m = R.modulo
    q,r = divmod(a,m)
    return QuotientRingElem{PolyRing{T},PolyRingElem{S,T}}(r, R)
end

###### Printing ######
function show(io::IO, R::QuotientRing{T,S}) where {T <: CRing, S <: CRingElem}
    print(io, "$(R.parent)/<$(R.modulo)>")
end

function show(io::IO, a::QuotientRingElem{T,S}) where {T <: CRing, S <: CRingElem}
    print(io, "$(a.d)")
end

###### Helper functions #######

function zero(R::QuotientRing{T,S}) where {T <: CRing, S <: CRingElem}
    return R(zero(R.parent))
end

function one(R::QuotientRing{T,S}) where {T <: CRing, S <: CRingElem}
    return R(one(R.parent))
end

function iszero(a::QuotientRingElem{T,S}) where {T <: CRing, S <: CRingElem}
    return iszero(a.d)
end

function isone(a::QuotientRingElem{T,S}) where {T <: CRing, S <: CRingElem}
    return isone(a.d)
end

function order(R::QuotientRing{PolyRing{T},PolyRingElem{S,T}}) where {T <: CRing, S <: CRingElem}
    return order(R.parent.parent)^degree(R.modulo)
end

###### Basic operators #######
function -(a::QuotientRingElem{T,S}) where {T <: CRing, S <: CRingElem}
    if iszero(a)
        return a
    else
        R = a.ring
        return R(-a.d)
    end
end

function +(a::QuotientRingElem{T,S}, b::QuotientRingElem{T,S}) where {T <: CRing, S <: CRingElem}
    if !sameRing(a, b)
        error("Cannot add elements of different rings")
    end
    z = a.d + b.d
    R = a.ring
    return R(z)
end

function -(a::QuotientRingElem{T,S}, b::QuotientRingElem{T,S}) where {T <: CRing, S <: CRingElem}
    return a+(-b)
end

function *(a::QuotientRingElem{T,S}, b::QuotientRingElem{T,S}) where {T <: CRing, S <: CRingElem}
    if !sameRing(a, b)
        error("Cannot multiply elements of different rings")
    end
    z = a.d * b.d
    R = a.ring
    return R(z)
end

function *(n::Integer, a::QuotientRingElem{T,S}) where {T <: CRing, S <: CRingElem}
    z = n*a.d
    R = a.ring
    return R(z)
end

function *(a::QuotientRingElem{T,S}, n::Integer) where {T <: CRing, S <: CRingElem}
    return n*a
end

function inv(x::QuotientRingElem{T,S}) where {T <: CRing, S <: CRingElem}
    R = x.ring
    g, a, b = gcdx(x.d, R.modulo)
    (!isone(g)) && error("Element $x is not invertible in $(x.ring)")
    return R(a)
end

#yet again exponensiation uses the simple square and multiply algorithm
function ^(x::QuotientRingElem{T,S}, y::U) where {T <: CRing, S <: CRingElem, U <: Integer}
    if iszero(x) || y == 1
        return x
    end
    R = x.ring
    if y == 0
        return one(R)
    end
    if y < 0
        x = inv(x)
        y = -y
    end
    # Should include some reduction of y modulo order if order is known
    bit = U(1) << (ndigits(y, base = 2) - 1)
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

function ==(a::QuotientRingElem{T,S}, b::QuotientRingElem{T,S}) where {T <: CRing, S <: CRingElem}
    return sameRing(a, b) && (a.d == b.d)
end


##########################################
#
#      Finite field of order p^2
#
##########################################
# Just a special case of the construction above for easy usage
# If p = 3 (mod 4), then -1 is a quadratic non-residue and x^2 + 1 is irredusible
function Fpsquare(p::T) where T <: Integer
    Zp = Zmod(p)
    ZpPol = PolyRing(Zp, "i")
    if !isprime(p)
        error("$p is not a prime!")
    end
    if p%4 != 3
        error("Only supported for p ≡ 3 (mod 4)")
    end
    m = ZpPol([Zp(T(1)),Zp(T(0)),Zp(T(1))])
    return QuotientRing(ZpPol, m)
end

# Element-Constructor
function (R::QuotientRing{PolyRing{Zmod{T}},PolyRingElem{ZmodElem{T},Zmod{T}}})(a::Array{S}) where {T <: Integer, S <: Integer}
    a = [T(x) for x in a]
    a = R.parent([R.parent.parent(i) for i in a])
    m = R.modulo
    q,r = divmod(a,m)
    return QuotientRingElem{PolyRing{Zmod{T}},PolyRingElem{ZmodElem{T},Zmod{T}}}(r, R)
end

#################################
#
# Dangerous or specific functions
#
#################################

function squareroot(a::ZmodElem{T}) where {T <: Integer}
    return a.ring(TShanks(a.d, a.ring.n)[1])
end

# This is just for the special case, Zp[X]/<f> where f is like x^2+1...
function cuberootOfOne(R::QuotientRing{PolyRing{Zmod{T}},PolyRingElem{ZmodElem{T},Zmod{T}}}) where {T <: Integer, S <: Integer}
    #idea is that the non trivial cube roots of 1 are -1±sqrt(-3)/2 as x^3-1 = (x-1)(x^2+x+1)
    Z = R.parent.parent
    if isone(Z(-3)^((Z.n-1)÷2))
        d = squareroot(Z(-3))
        u = R([d])
    elseif isone(Z(3)^((Z.n-1)÷2))
        d = squareroot(Z(3))
        u = R([zero(Z), d])
    else
        error("Probably no cuberoot of 1 in the field")
    end
    id = one(R)
    z = inv(2*id)
    return (-id+u)*z
end

# Make sure that the ring is actually a field before using this... Further 3|(q-1)
function cuberoot(a::CRingElem)
    q = order(a.ring)
    if q%9 != 7
        error("Currently only supported for q ≡ 7 (mod 9)")
    end
    #general euler criterion
    k = (q-1)÷3
    if !isone(a^k)
        error("$a is not a cube in $(a.ring)")
    end
    k = (q+2)÷9
    return a^k
end

function allcuberoots(a::CRingElem)
    w = cuberootOfOne(a.ring)
    k = cuberoot(a)
    return (k, k*w, k*w^2)
end

function allElements(R::QuotientRing{PolyRing{Zmod{T}},PolyRingElem{ZmodElem{T},Zmod{T}}}) where T <: Integer
    deg = degree(R.modulo)
    if deg != 2
        println("Too many!")
    end
    Z = R.parent.parent
    l = QuotientRingElem{PolyRing{Zmod{T}},PolyRingElem{ZmodElem{T},Zmod{T}}}[]
    for i1 in allElements(Z)
        for i2 in allElements(Z)
            push!(l,R([i1,i2]))
        end
    end
    return l
end


#################################
#
#            Other
#
#################################

legendre(a, p) = powermod(a, (p - 1) ÷ 2, p)

function TShanks(n::T, p::T) where T <: Union{Int, Int128, BigInt}
    legendre(n, p) != 1 && throw(ArgumentError("$n not a square (mod $p)"))
    local q::T = p - one(p)
    local s::T = 0
    while iszero(q % 2)
        q ÷= 2
        s += one(s)
    end
    if s == one(s)
        r = powermod(n, (p + 1) >> 2, p)
        return r, p - r
    end
    local z::T
    for z in 2:(p - 1)
        p - 1 == legendre(z, p) && break
    end
    local c::T = powermod(z, q, p)
    local r::T = powermod(n, (q + 1) >> 1, p)
    local t::T = powermod(n, q, p)
    local m::T = s
    local t2::T = zero(p)
    while !iszero((t - 1) % p)
        t2 = (t * t) % p
        local i::T
        for i in Base.OneTo(m)
            iszero((t2 - 1) % p) && break
            t2 = (t2 * t2) % p
        end
        b = powermod(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    end
    return r, p - r
end

end
