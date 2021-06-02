module TwistedHessianEC

export EllipticCurve, EllipticCurveElem, TwiHesEC, TwiHesECElem, identity, isIdentity, normalized,
jInvariant, Isogeny, allPointsNaive

include("FieldArithmetic.jl")

using .FieldArithmetic
import Base: +,-,*,/,^,==,inv,show,zeros,gcd,gcdx,iszero,isone,zero,one,identity
import .FieldArithmetic: order

export Ring, RingElem, CRing, CRingElem, Field, FieldElem, Zmod, ZmodElem, PolyRing, PolyRingElem,
QuotientRing, QuotientRingElem, Fpsquare, order, squareroot, cuberoot, cuberootOfOne, allcuberoots,
allElements

abstract type EllipticCurve end
abstract type EllipticCurveElem end

struct TwiHesEC{T <: CRing, S <: CRingElem} <: EllipticCurve
    parent::T
    a::S
    d::S
end

struct TwiHesECElem{T <: CRing, S <: CRingElem} <: EllipticCurveElem
    X::S
    Y::S
    Z::S
    curve::TwiHesEC{T,S}
end

# Curve-Constructor
function TwiHesEC(R::T, a::S, d::S) where {T <: CRing, S <: CRingElem}
    if isSingular(R,a,d)
        error("Curve is singular!")
    end
    return TwiHesEC{T,S}(R, a, d)
end

# Element-Constructor
function (E::TwiHesEC{T,S})(X::S, Y::S, Z::S) where {T <: CRing, S <: CRingElem}
    if !validPoint(X, Y, Z, E)
        error("Point ([$(X)] : [$(Y)] : [$(Z)]) is not valid!!!")
    end
    return TwiHesECElem{T,S}(X, Y, Z, E)
end

###### Printing ######
function show(io::IO, E::TwiHesEC{T,S}) where {T <: CRing, S <: CRingElem}
    print(io, "Elliptic curve [$(E.a)]*X^3 + Y^3 + Z^3 = [$(E.d)]*X*Y*Z defined over $(E.parent)")
end

function show(io::IO, P::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    print(io, "([$(P.X)] : [$(P.Y)] : [$(P.Z)])")
end

###### Helper-functions ######
function isSingular(E::T, a::S, d::S) where {T <: CRing, S <: CRingElem}
    return (iszero(a*(27*a - d^3)))
end

function validPoint(X::S, Y::S, Z::S, E::TwiHesEC{T,S}) where {T <: CRing, S <: CRingElem}
    return E.a*X^3 + Y^3 + Z^3 == E.d*X*Y*Z && !isZeroCoordinate(X,Y,Z)
end

function isZeroCoordinate(X::S, Y::S, Z::S) where {S <: CRingElem}
    F = X.ring
    return (X,Y,Z)==(zero(F),zero(F),zero(F))
end

function identity(E::TwiHesEC{T,S}) where {T <: CRing, S <: CRingElem}
    F = E.parent
    return E(zero(F),-one(F),one(F))
end

function isIdentity(P::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    return P == identity(P.curve)
end

function sameCurve(P::TwiHesECElem{T,S}, Q::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    return P.curve == Q.curve
end

function normalized(P::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    X, Y, Z = P.X, P.Y, P.Z
    E = P.curve
    F = E.parent
    if !iszero(Z)
        Zinv = inv(Z)
        return E(X*Zinv, Y*Zinv, one(F))
    else
        Xinv = inv(X)
        return E(one(F), Y*Xinv, zero(F))
    end
end

#Snail speed check...
function order(P::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    i = 1
    Pi = P
    while !isIdentity(Pi)
        Pi = Pi + P
        i += 1
    end
    return i
end

###### Basic operators #######

function -(P::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    E = P.curve
    return E(P.X, P.Z, P.Y)
end

function +(P::TwiHesECElem{T,S}, Q::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    if !(sameCurve(P,Q))
        error("Cannot add points from different curves!")
    end
    E = P.curve
    X1, Y1, Z1 = P.X, P.Y, P.Z
    X2, Y2, Z2 = Q.X, Q.Y, Q.Z
    X = X1^2*Y2*Z2 - X2^2*Y1*Z1
    Y = Z1^2*X2*Y2 - Z2^2*X1*Y1
    Z = Y1^2*X2*Z2 - Y2^2*X1*Z1
    if isZeroCoordinate(X,Y,Z)
        return fallbackAdd(P,Q)
    else
        return E(X, Y, Z)
    end
end

function fallbackAdd(P::TwiHesECElem{T,S}, Q::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    E = P.curve
    a = E.a
    X1, Y1, Z1 = P.X, P.Y, P.Z
    X2, Y2, Z2 = Q.X, Q.Y, Q.Z
    X = Z2^2*X1*Z1 - Y1^2*X2*Y2
    Y = Y2^2*Y1*Z1 - a*X1^2*X2*Z2
    Z = a*X2^2*X1*Y1 - Z1^2*Y2*Z2
    return E(X, Y, Z)
end

function -(P::TwiHesECElem{T,S}, Q::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    return P+(-Q)
end

#Regard E as a Z-module... This is again square and multiply (or strictly speaking double and add)
function *(n::U, P::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem, U <: Integer}
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
        Q = Q+Q
        if (bit & n) != 0
            Q += P
        end
        bit >>=1
    end
    return Q
end

function *(P::TwiHesECElem{T,S}, n::U) where {T <: CRing, S <: CRingElem, U <: Integer}
    return n*P
end

function ==(P::TwiHesECElem{T,S}, Q::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    P1 = normalized(P)
    Q1 = normalized(Q)
    return sameCurve(P,Q) && P1.X == Q1.X && P1.Y == Q1.Y && P1.Z == Q1.Z
end

############################
#
#   Elliptic Curve Stuff
#
############################

function jInvariant(E::TwiHesEC{T,S}) where {T <: CRing, S <: CRingElem}
    a = E.a
    d = E.d
    over = -(216*a*d+d^4)^3
    under = a*(27*a-d^3)^3
    return over*inv(under)
end

############################
#
#        Isogenies
#
############################

struct Isogeny
    domain::EllipticCurve
    codomain::EllipticCurve
    carryP::Function
    carryPfallback::Function
end

# Constructor
function Isogeny(E::EllipticCurve, P::EllipticCurveElem)
    if isIdentity(P)
        error("$P is the identity!")
    end
    if isIdentity(2*P)
        return deg2_isogeny(E,P)
    end
    if isIdentity(3*P)
        return deg3_isogeny(E,P)
    end
    error("$P is not of order 2 or 3")
end

#Later optimization: Normalize is only needed for computing the D of the isogenous
#curve. This is not necessary until the very last step. Maybe we can skip it?
#NB! Then formulas change, so be a bit careful (this one relies on Z = 1).
#Reasoning for the last part: Thm 1 from Dang and Moody says that order(P) = 3 <=> at one coord of P is zero
function deg2_isogeny(E::TwiHesEC{T,S}, P::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    P1 = normalized(P)
    a = E.a
    d = E.d
    F = E.parent
    A = a^2
    s = P1.X
    t = P1.Y
    D = (-d+6*inv(s*t))*inv(s)
    E1 = TwiHesEC(F, A, D)
    phi(X, Y, Z) = (X*(X^2*t - s^2*Y*Z), Y*(Z^2*s*t - X*Y), Z*(Y^2*s - t^2*X*Z))
    phi2(X, Y, Z) =(X*(X*Z - Y^2*s*t), Y*(Y*Z*t^2 - X^2*a*s), Z*(X*Y*a*s^2 - Z^2*t))
    return Isogeny(E, E1, phi, phi2)
end

#Two different isogenies based on which subgroup we are looking at
function deg3_isogeny(E::TwiHesEC{T,S}, P::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    if iszero(P.X)
        return deg3_isogeny_c1(E, P)
    else
        return deg3_isogeny_c2(E, P)
    end
end

#Case 1: {(0:-1:1),(0:-w:1),(0:-w^2:1)}, with w^3 = 1, w != 1
function deg3_isogeny_c1(E::TwiHesEC{T,S}, P::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    a = E.a
    d = E.d
    F = E.parent
    A = d^3-27*a
    D = 3*d
    E1 = TwiHesEC(F, A, D)
    w = -normalized(P).Y
    phi(X, Y, Z) = (X*Y*Z, a*X^3 + w^2*Y^3 + w*Z^3, a*X^3 + w*Y^3 + w^2*Z^3)
    return Isogeny(E, E1, phi, phi)
end

#Case 2: {(0:-1:1),(1:-c:0),(1:0:-c)}, with c^3 = a
function deg3_isogeny_c2(E::TwiHesEC{T,S}, P::TwiHesECElem{T,S}) where {T <: CRing, S <: CRingElem}
    P1 = normalized(P)
    a = E.a
    d = E.d
    F = E.parent
    if iszero(P1.Z)
        c = -P1.Y
    else
        c = inv(-P1.X)
    end
    A = d^2*c + 3*d*c^2 + 9*a
    D = d + 6*c
    E1 = TwiHesEC(F, A, D)
    phi(X, Y, Z) = (X*Y*Z, c^2*X^2*Z + c*X*Y^2 + Y*Z^2, c^2*X^2*Y + c*X*Z^2 + Y^2*Z)
    return Isogeny(E, E1, phi, phi)
end



###### Printing ######
function show(io::IO, phi::Isogeny)
    print(io, "Isogeny φ : E_1 -> E_2 where \nE_1 = $(phi.domain)\nE_2 = $(phi.codomain)")
end

###### Go brrr ######
# Evaluate a point
function (phi::Isogeny)(P::EllipticCurveElem)
    if P.curve != phi.domain
        error("Point $P is not from the right curve")
    end
    (X, Y, Z) = phi.carryP(P.X,P.Y,P.Z)
    if !isZeroCoordinate(X, Y, Z)
        return phi.codomain(X, Y, Z)
    else
        (X, Y, Z) = phi.carryPfallback(P.X,P.Y,P.Z)
        return phi.codomain(X, Y, Z)
    end
end


############################
#
#           Other
#
############################
function allPointsNaive(E::TwiHesEC{T,S}) where {T <: CRing, S <: CRingElem}
    points = TwiHesECElem{T,S}[]
    Z = one(E.parent)
    for X in allElements(E.parent)
        for Y in allElements(E.parent)
            if validPoint(X, Y, Z, E)
                push!(points,E(X, Y, Z))
            end
        end
    end
    Z = zero(E.parent)
    X = one(E.parent)
    for Y in allElements(E.parent)
        if validPoint(X, Y, Z, E)
            push!(points,E(X, Y, Z))
        end
    end
    return points
end

end
