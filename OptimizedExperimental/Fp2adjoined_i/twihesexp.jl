include("fpsquarei.jl")

# A projective point on the curve.
# The curve itself is not explicitly defined
struct TwiHesPoint{T <: Integer}
    X::Fp2Elem{T}
    Y::Fp2Elem{T}
    Z::Fp2Elem{T}
end

# Element-Constructor
# NB! Assumes point is valid! Checking everytime you construct a new point is very slow....
function EC(X::Fp2Elem{T}, Y::Fp2Elem{T}, Z::Fp2Elem{T}) where {T <: Integer}
    return TwiHesPoint{T}(X, Y, Z)
end

############################
#
# Printing
#
############################

function show(io::IO, P::TwiHesPoint{T}) where {T <: Integer}
    print(io, "($(P.X) : $(P.Y) : $(P.Z))")
end

############################
#
# Some helper functions
#
############################

function isValid(P::TwiHesPoint{T}, a::Fp2Elem{T}, d::Fp2Elem{T}) where {T <: Integer}
    return a*P.X^3 + P.Y^3 + P.Z^3 == d*P.X*P.Y*P.Z && !isZeroCoordinate(P.X,P.Y,P.Z)
end

function isZeroCoordinate(X::Fp2Elem{T}, Y::Fp2Elem{T}, Z::Fp2Elem{T}) where {T <: Integer}
    return (X,Y,Z)==(zero(X.field),zero(X.field),zero(X.field))
end

function identity(F::Fp2{T}) where {T <: Integer}
    return EC(zero(F),-one(F),one(F))
end

function isIdentity(P::TwiHesPoint{T}) where {T <: Integer}
    return P == identity(P.X.field)
end

function ==(P::TwiHesPoint{T}, Q::TwiHesPoint{T}) where {T <: Integer}
    P1 = normalized(P)
    Q1 = normalized(Q)
    return P1.X == Q1.X && P1.Y == Q1.Y && P1.Z == Q1.Z
end

# Prefers Z = 1, then Y = 1
function normalized(P::TwiHesPoint{T}) where {T <: Integer}
    if !iszero(P.Z)
        Zinv = inv(P.Z)
        return EC(P.X*Zinv, P.Y*Zinv, one(P.X.field))
    else
        return EC(P.X*inv(P.Y), one(P.X.field), zero(P.X.field))
    end
end

# Recover curve coefficient d from a and a point P, with [3]P != O
function recover_d(P::TwiHesPoint{T}, a::Fp2Elem{T}) where {T <: Integer}
    return (a*P.X*P.X*P.X + P.Y*P.Y*P.Y + P.Z*P.Z*P.Z)*inv(P.X*P.Y*P.Z)
end

############################
#
# Basic operators
#
############################

# Inverse point
# Input: A point P
# Output: -P
function -(P::TwiHesPoint{T}) where {T <: Integer}
    return EC(P.X, P.Z, P.Y)
end

# Addition based on the rotated addition law.
# Input: Points P, Q on curve E(a,d), and curve-parameter a
# Output: P+Q
function add(P::TwiHesPoint{T}, Q::TwiHesPoint{T}, a::Fp2Elem{T}) where {T <: Integer}
    A = P.X*Q.Z
    B = P.Z*Q.Z
    C = P.Y*Q.X
    D = P.Y*Q.Y
    E = P.Z*Q.Y
    F = a*P.X*Q.X
    G = (D + B)*(A - C)
    H = (D - B)*(A + C)
    J = (D + F)*(A - E)
    K = (D - F)*(A + E)
    X = G - H
    Y = K - J
    Z = J + K - G - H - 2*(B - F)*(C + E)
    return EC(X, Y, Z)
end

# This function uses the simple double and add.
# Input: integer n and point P
# Output: [n]P
function mul(n::Integer, P::TwiHesPoint{T}, a::Fp2Elem{T}) where {T <: Integer}
    if n == 1
        return P
    end
    if n == 0
        return identity(a.field)
    end
    if n < 0
        P = -P
        n = -n
    end
    bit = T(1) << (ndigits(n, base = 2) - 1)
    Q = P
    bit >>= 1
    while bit != 0
        Q = double(Q, a)
        if (bit & n) != 0
            Q = add(Q, P, a)
        end
        bit >>=1
    end
    return Q
end

# Doubling of a point
# Costs 6M + 3S + 1Ma
# Can be improved, but we wish to avoid the need for d
# Input: A point P on E(a,d), and curve parameter a
# Output: [2]P
function double(P::TwiHesPoint{T}, a::Fp2Elem{T}) where {T <: Integer}
    X, Y, Z = P.X, P.Y, P.Z
    #X = Z^3*X - Y^3*X
    #Y = Y^3*Z - a*X^3*Z
    #Z = a*X^3*Y - Z^3*Y
    X3 = P.X*P.X*P.X
    Y3 = P.Y*P.Y*P.Y
    Z3 = P.Z*P.Z*P.Z
    aX3 = a*X3
    X_new = (Z3 - Y3)*P.X
    Y_new = (Y3 - aX3)*P.Z
    Z_new = (aX3 - Z3)*P.Y
    return EC(X_new, Y_new, Z_new)
end

#Repeated doubling
#Input point P, integer e
#Output [2^e]P
function doublePow(P::TwiHesPoint{T}, a::Fp2Elem{T}, e::Integer) where {T <: Integer}
    Q = P
    for i in 1:e
        Q = double(Q, a)
    end
    return Q
end

############################
#
# Isogeny-related
#
############################

#Input: Curve a and d from E(a,d) in twisted Hessian form
#Output: j(E)
function jInvariant(a::Fp2Elem{T}, d::Fp2Elem{T}) where {T <: Integer}
    d3 = d*d*d
    o = (216*a+d3)*d
    u = (d3-27*a)
    over = o*o*o
    under = a*u*u*u
    return over*inv(under)
end

#input: A point P, and precomputed values a2 = a^2, b2 = b^2 ab = a*b from the kernel (a,b,b)
#Output: Phi(P)
function iso2eval(P::TwiHesPoint{T}, a2::Fp2Elem{T}, b2::Fp2Elem{T}, ab::Fp2Elem{T}) where {T <: Integer}
    X2 = P.X*P.X
    Y2 = P.Y*P.Y
    Z2 = P.Z*P.Z
    XY = P.X*P.Y
    XZ = P.X*P.Z
    YZ = P.Y*P.Z
    X_new = P.X*(X2*b2 - a2*YZ)
    Y_new = P.Y*(Z2*ab - b2*XY)
    Z_new = P.Z*(Y2*ab - b2*XZ)
    return EC(X_new, Y_new, Z_new)
end

# Some precomputed values for faster 2-isogenies
# Input: Point G of order 2
# Output: Precomputed values for evaluating isogeny generated by G
function iso2precomp(G::TwiHesPoint{T}) where {T <: Integer}
    return G.X*G.X, G.Y*G.Y, G.X*G.Y
end

# 3-isogeny evaluation when the kernel is <1 : -c : 0)>
# Input: Point P, and cb and b such that c^3 = a
# Output: Phi(P)
function iso3eval_c2(P::TwiHesPoint{T}, c::Fp2Elem{T}) where {T <: Integer}
    YZ = P.Y*P.Z
    cX = c*P.X
    c2X2 = cX*cX
    X_new = P.X*YZ
    Y_new = c2X2*P.Z + cX*P.Y*P.Y + YZ*P.Z
    Z_new = c2X2*P.Y + cX*P.Z*P.Z + YZ*P.Y
    return EC(X_new, Y_new, Z_new)
end

# Evaluating the isomorphism \varphi: H(a, d) -> H(1, d/a^(1/3))
# Input A point P, and a cuberoot of A
# Output \varphi(P)
function isoHessianEval(P::TwiHesPoint{T}, cuberoota::Fp2Elem{T}) where {T <: Integer}
    X, Y, Z = P.X, P.Y, P.Z
    return EC(cuberoota*X, Y, Z)
end

###############################
#
# Optimal isogeny strategies
#
###############################

#Input: A generator point G of order 2^e, and two other points P, Q, curve parameter a, isogeny degree e, and an optimal strategy S
#Output: E/<G>, and phi(P), phi(Q)
function iso2pow_opt(G::TwiHesPoint{T}, P::TwiHesPoint{T}, Q::TwiHesPoint{T}, a::Fp2Elem{T}, e::Integer, Strat::Vector{Int64}) where {T <: Integer}
    Pi = P
    Qi = Q
    ai = a
    Si = Vector{Tuple{Integer,TwiHesPoint{T}}}()
    push!(Si, (e, G))
    i = 1
    while !isempty(Si)
        (h, Ki) = pop!(Si)
        if h == 1
            ai = ai*ai
            a2, b2, ab = iso2precomp(Ki)
            S_new = Vector{Tuple{Integer,TwiHesPoint{T}}}()
            while !isempty(Si)
                (h, Xi) = popfirst!(Si)
                push!(S_new, (h-1, iso2eval(Xi, a2, b2, ab)))
            end
            Si = S_new
            Pi = iso2eval(Pi, a2, b2, ab)
            Qi = iso2eval(Qi, a2, b2, ab)
        else
            strat_i = Strat[i]
            if 0 < strat_i < h
                push!(Si, (h, Ki))
                push!(Si, (h-strat_i, doublePow(Ki, ai, strat_i)))
                i += 1
            else
                error("Strategy seems to be invalid")
            end
        end
    end
    return Pi, Qi, ai
end

#Second part, no points to carry
#Input: A generator point G of order 2^e, curve parameter a, isogeny degree e, and an optimal strategy S
#Output: parameters a' and d' of ending curve
function iso2pow_opt(G::TwiHesPoint{T}, a::Fp2Elem{T}, e::Integer, Strat::Vector{Int64}) where {T <: Integer}
    ai = a
    Si = Vector{Tuple{Integer,TwiHesPoint{T}}}()
    push!(Si, (e, G))
    i = 1
    while !isempty(Si)
        (h, Ki) = pop!(Si)
        if h == 1
            if isempty(Si)
                #Second to last curve, must recover d here, and use it to find d on the final curve
                di = recover_d(Ki, ai)
                a_final = ai*ai
                d_final = Ki.Y*(6*Ki.Y-Ki.X*di)*inv(Ki.X*Ki.X)
                return a_final, d_final
            end
            ai = ai*ai
            a2, b2, ab = iso2precomp(Ki)
            S_new = Vector{Tuple{Integer,TwiHesPoint{T}}}()
            while !isempty(Si)
                (h, Xi) = popfirst!(Si)
                push!(S_new, (h-1, iso2eval(Xi, a2, b2, ab)))
            end
            Si = S_new
        else
            strat_i = Strat[i]
            if 0 < strat_i < h
                push!(Si, (h, Ki))
                push!(Si, (h-strat_i, doublePow(Ki, ai, strat_i)))
                i += 1
            else
                error("Strategy seems to be invalid")
            end
        end
    end
end


function Xnormalized(P::TwiHesPoint{T}) where {T <: Integer}
    X, Y, Z = P.X, P.Y, P.Z
    F = X.field
    if !iszero(X)
        Xinv = inv(X)
        return EC(one(F), Y*Xinv, Z*Xinv)
    else
        Zinv = inv(Z)
        return EC(zero(F), Y*Zinv, one(F))
    end
end