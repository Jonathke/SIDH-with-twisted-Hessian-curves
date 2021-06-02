##########################################
#         Importing
##########################################

include("FieldArithmetic.jl")

using .FieldArithmetic
import Base: +,-,*,/,^,==,inv,show,zeros,gcd,gcdx,iszero,isone,zero,one,identity
import .FieldArithmetic: order

abstract type EllipticCurve end
abstract type EllipticCurveElem end

##########################################
#
#Elliptic Curve in twisted Hessian form
#
##########################################

struct Affine_EC_TwiHes{T <: CRing, S <: CRingElem} <: EllipticCurve
    parent::T
    a::S
    d::S
end

struct Affine_EC_TwiHesElem{T <: CRingElem, S <: CRing} <: EllipticCurveElem
    x::T
    y::T
    curve::Affine_EC_TwiHes{S, T}
end

# Curve-Constructor
function Affine_EC_TwiHes(R::T, a::S, d::S) where {T <: CRing, S <: CRingElem}
    if isSingular(R,a,d)
        error("Curve is singular!")
    end
    return Affine_EC_TwiHes{T,S}(R, a, d)
end

# Element-Constructor
function (E::Affine_EC_TwiHes{T,S})(x::S, y::S) where {T <: CRing, S <: CRingElem}
    if !validPoint(x, y, E)
        error("Point ($x, $y) is not valid!!!")
    end
    return Affine_EC_TwiHesElem{S, T}(x,y,E)
end


###### Printing ######
function show(io::IO, R::Affine_EC_TwiHes{T,S}) where {T <: CRing, S <: CRingElem}
    print(io, "Elliptic curve [$(R.a)]*x^3 + y^3 + 1 = [$(R.d)]*x*y defined over $(R.parent)")
end

function show(io::IO, P::Affine_EC_TwiHesElem{T,S}) where {T <: CRingElem, S <: CRing}
    print(io, "([$(P.x)], [$(P.y)])")
end

###### Helper-functions ######
function isSingular(R::T, a::S, d::S) where {T <: CRing, S <: CRingElem}
    return (iszero(a*(27*a - d^3)))
end

function validPoint(x::S, y::S, E::Affine_EC_TwiHes{T,S}) where {T <: CRing, S <: CRingElem}
    return E.a*x^3 + y^3 + one(E.parent) == E.d*x*y
end

function identity(E::Affine_EC_TwiHes{T,S}) where {T <: CRing, S <: CRingElem}
    F = E.parent
    return E(zero(F),-one(F))
end

function isIdentity(P::Affine_EC_TwiHesElem{T,S}) where {T <: CRingElem, S <: CRing}
    return P == identity(P.curve)
end

function sameCurve(P::Affine_EC_TwiHesElem{T,S}, Q::Affine_EC_TwiHesElem{T,S}) where {T <: CRingElem, S <: CRing}
    return (P.curve == Q.curve)
end

#Snail speed, but oh well
function order(P::Affine_EC_TwiHesElem{T,S}) where {T <: CRingElem, S <: CRing}
    i = 1
    Pi = P
    while !isIdentity(Pi)
        Pi = Pi + P
        i += 1
    end
    return i
end

###### Basic operators #######

#How to invert (x,0)??
#Answer: A point (x,0) implies that a*x^3 = -1, or x = a^(-1/3), which is only possible if a is a cube.
function -(P::Affine_EC_TwiHesElem{T,S}) where {T <: CRingElem, S <: CRing}
    E = P.curve
    return E(P.x*inv(P.y),inv(P.y))
end

function +(P::Affine_EC_TwiHesElem{T,S}, Q::Affine_EC_TwiHesElem{T,S}) where {T <: CRingElem, S <: CRing}
    if !(sameCurve(P, Q))
        error("Cannot add points from different curves!")
    end
    if P == Q
        return double(P)
    end
    E = P.curve
    a = E.a
    x1, y1, x2, y2 = P.x, P.y, Q.x, Q.y
    z1 = x1 - y1^2 * x2 * y2
    z2 = y1 * y2^2 - a * x1^2 * x2
    z3 = inv(a * x1 * y1 * x2^2 - y2)
    return E(z1*z3, z2*z3)
end

function -(P::Affine_EC_TwiHesElem{T,S}, Q::Affine_EC_TwiHesElem{T,S}) where {T <: CRingElem, S <: CRing}
    return P+(-Q)
end

function double(P::Affine_EC_TwiHesElem{T,S}) where {T <: CRingElem, S <: CRing}
    E = P.curve
    a = E.a
    x, y = P.x, P.y
    z1 = x - y^3 * x
    z2 = y^3 - a * x^3
    z3 = inv(a * y * x^3 - y)
    return E(z1*z3, z2*z3)
end

#Regard E as a Z-module... This is again square and multiply (or strictly speaking double and add)
function *(n::U, P::Affine_EC_TwiHesElem{T,S}) where {T <: CRingElem, S <: CRing, U <: Integer}
    if isIdentity(P) || n == 1
        return P
    end
    E = P.curve
    if n == 0
        return identity(E)
    end
    if n < 0
        P = -P
        n = -n
    end
    bit = U(1) << (ndigits(n, base = 2) - 1)
    Q = P
    bit >>= 1
    while bit != 0
        Q = double(Q)
        if (bit & n) != 0
            Q += P
        end
        bit >>=1
    end
    return Q
end

function *(P::Affine_EC_TwiHesElem{T,S}, n::U) where {T <: CRingElem, S <: CRing, U <: Integer}
    return n*P
end

function ==(P::Affine_EC_TwiHesElem{T,S}, Q::Affine_EC_TwiHesElem{T,S}) where {T <: CRingElem, S <: CRing}
    return sameCurve(P,Q) && (P.x == Q.x) && (P.y == Q.y)
end

############################
#
#   Elliptic Curve Stuff
#
############################

#See notebook for calculation. It is calculated as the j-invariant of the isomorphic Weierstrass-curve
function jInvariant(E::Affine_EC_TwiHes{T,S}) where {T <: CRing, S <: CRingElem}
    a = E.a
    d = E.d
    over = -(216*a*d+d^4)^3
    under = a*(27*a-d^3)^3
    return over*inv(under)
end

function toWeierstrass(E::Affine_EC_TwiHes{T,S}) where {T <: CRing, S <: CRingElem}
    a = E.a
    d = E.d
    w1 = -8*a*d - d^4*inv(27*one(E.parent))
    w2 = 16*a^2 + 40*a*d^3*inv(27*one(E.parent)) - 2*d^6*inv(729*one(E.parent))
    println(E)
    println("is isomorphic (over the algebraic clousure) to the curve")
    println("y^2 = x^3 + [$w1]*x + [$w2]")
end

#FromWeierstrass: Find point of order 3, then transform to triangular form and finally to twisted Hessian

############################
#
#           Other
#
############################
function allPointsNaive(E::Affine_EC_TwiHes{T,S}) where {T <: CRing, S <: CRingElem}
    points = Affine_EC_TwiHesElem{S,T}[]
    for x in allElements(E.parent)
        for y in allElements(E.parent)
            if validPoint(x, y, E)
                push!(points,E(x,y))
            end
        end
    end
    return points
end

#Disregard
function WeierjInvariant(a, b)
    delta = -16*(4*a^3+27*b^2)
    j = -1728*(4*a)^3*inv(delta)
    println(j)
end


############################
#
#        Isogenies
#
############################

struct Affine_isogeny
    domain::EllipticCurve
    codomain::EllipticCurve
    carryP::Function
end

# Constructor
function Affine_isogeny(E::EllipticCurve, P::EllipticCurveElem)
    if isIdentity(2*P)
        return deg2_isogeny(E,P)
    end
    if isIdentity(3*P)
        return deg3_isogeny(E,P)
    end
    error("$P is not of order 2 or 3")
end

function deg3_isogeny(E::EllipticCurve, P::EllipticCurveElem)

end

function deg2_isogeny(E::Affine_EC_TwiHes, P::Affine_EC_TwiHesElem)
    a = E.a
    d = E.d
    F = E.parent
    A = a^2
    s = P.x
    t = P.y
    D = (-d+6*inv(s*t))*inv(s)
    E1 = Affine_EC_TwiHes(F, A, D)
    phi(x,y) = (x*((x-s*t*y^2)*inv(a*s^2*x*y-t)),y*((y*t^2-a*s*x^2)*inv(a*s^2*x*y-t)))
    return Affine_isogeny(E, E1, phi)
end


# Evaluate a point
function (phi::Affine_isogeny)(P::EllipticCurveElem)
    if P.curve != phi.domain
        error("Point $P is not from the right curve")
    end
    (x,y) = phi.carryP(P.x,P.y)
    return phi.codomain(x,y)
end


###### Printing ######
function show(io::IO, phi::Affine_isogeny)
    print(io, "Isogeny φ : E_1 -> E_2 where \nE_1 = $(phi.domain)\nE_2 = $(phi.codomain)")
end
