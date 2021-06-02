##########################################
#         Importing
##########################################

include("FieldArithmetic.jl")

using .FieldArithmetic
import Base: +,-,*,/,^,==,inv,show,zeros,gcd,gcdx,iszero,isone,zero,one

##########################################
#
#Elliptic Curve in Weierstrass Normal Form
#
##########################################

struct Weierstrass_EC{T <: CRing, S <: CRingElem}
    parent::T
    A::S
    B::S
end

struct Weierstrass_ECElem{T <: CRingElem, S <: CRing}
    X::T
    Y::T
    Z::T
    curve::Weierstrass_EC{S, T}
end

# Curve-Constructor
function Weierstrass_EC(R::T, A::S, B::S) where {T <: CRing, S <: CRingElem}
    if isSingular(R,A,B)
        error("Curve is singular!")
    end
    return Weierstrass_EC{T,S}(R, A, B)
end

# Element-Constructor
function (E::Weierstrass_EC{T,S})(X::S, Y::S, Z::S) where {T <: CRing, S <: CRingElem}
    if !validPoint(X, Y, Z, E)
        error("Point is not valid!!!")
    end
    return Weierstrass_ECElem{S, T}(X,Y,Z,E)
end


###### Printing ######
function show(io::IO, R::Weierstrass_EC{T,S}) where {T <: CRing, S <: CRingElem}
    print(io, "Elliptic curve y^2 = x^3 + $(R.A)*x + $(R.B) defined over $(R.parent)")
end

function show(io::IO, P::Weierstrass_ECElem{T,S}) where {T <: CRingElem, S <: CRing}
    print(io, "($(P.X): $(P.Y): $(P.Z))")
end

###### Helper-functions ######
function isSingular(E::T, A::S, B::S) where {T <: CRing, S <: CRingElem}
    return iszero(4*A^3 + 27*B^2)
end

function validPoint(X::S, Y::S, Z::S, E::Weierstrass_EC{T,S}) where {T <: CRing, S <: CRingElem}
    return X^3 + X*Z^2*E.A + Z^3*E.B == Z*Y^2 && !isZeroCoordinate(X, Y, Z)
end

function isZeroCoordinate(X::S, Y::S, Z::S) where {S <: CRingElem}
    F = X.ring
    return (X,Y,Z)==(zero(F),zero(F),zero(F))
end

function identity(E::Weierstrass_EC{T,S}) where {T <: CRing, S <: CRingElem}
    F = E.parent
    return E(zero(F),one(F),zero(F))
end

function isIdentity(P::Weierstrass_ECElem{S, T}) where {T <: CRing, S <: CRingElem}
    return P == identity(P.curve)
end

function sameCurve(P::Weierstrass_ECElem{S, T}, Q::Weierstrass_ECElem{S, T}) where {T <: CRing, S <: CRingElem}
    return P.curve == Q.curve
end

function normalized(P::Weierstrass_ECElem{S, T}) where {T <: CRing, S <: CRingElem}
    X, Y, Z = P.X, P.Y, P.Z
    E = P.curve
    F = E.parent
    if !iszero(Z)
        Zinv = inv(Z)
        return E(X*Zinv, Y*Zinv, one(F))
    else
        return identity(E)
    end
end

###### Basic operators #######

function -(P::Weierstrass_ECElem{S, T}) where {T <: CRing, S <: CRingElem}
    E = P.curve
    return E(P.X, -P.Y, P.Z)
end

function +(P::Weierstrass_ECElem{S, T}, Q::Weierstrass_ECElem{S, T}) where {T <: CRing, S <: CRingElem}
    if !(sameCurve(P,Q))
        error("Cannot add points from different curves!")
    end
    if isIdentity(P)
        return Q
    end
    if isIdentity(Q)
        return P
    end
    E = P.curve
    P = normalized(P)
    Q = normalized(Q)
    A = E.A
    X1, Y1, Z1 = P.X, P.Y, P.Z
    X2, Y2, Z2 = Q.X, Q.Y, Q.Z
    if X1 == X2 && Y1 == -Y2
        return identity(E)
    else
        if P == Q
            l = (3*X1^2 + A)*inv(2*Y1)
        else
            l = (Y2 - Y1)*inv(X2 - X1)
        end
        X3 = l^2 - X1 - X2
        Y3 = l*(X1 - X3) - Y1
        return E(X3, Y3, one(F))
    end
end

function -(P::Weierstrass_ECElem{S, T}, Q::Weierstrass_ECElem{S, T}) where {T <: CRing, S <: CRingElem}
    return P+(-Q)
end

function *(n::U, P::Weierstrass_ECElem{S, T}) where {T <: CRing, S <: CRingElem, U <: Integer}
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

function *(P::Weierstrass_ECElem{S, T}, n::U) where {T <: CRing, S <: CRingElem, U <: Integer}
    return n*P
end

function ==(P::Weierstrass_ECElem{S, T}, Q::Weierstrass_ECElem{S, T}) where {T <: CRing, S <: CRingElem}
    P1 = normalized(P)
    Q1 = normalized(Q)
    return sameCurve(P,Q) && P1.X == Q1.X && P1.Y == Q1.Y && P1.Z == Q1.Z
end
