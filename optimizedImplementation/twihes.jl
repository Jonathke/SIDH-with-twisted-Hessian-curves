include("fpsquare.jl")

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

# Tripling of a point
# This requires d != 0. While not ideal, we avoid it by starting "far away" from E with j(E) = 0
# Costs 8M + 6S + 1Ma + 2Md
# Input: A point P on E(a,d), and curve parameters a and d
# Output: [3]P
function triple(P::TwiHesPoint{T}, a::Fp2Elem{T}, d::Fp2Elem{T}) where {T <: Integer}
    X, Y, Z = P.X, P.Y, P.Z
    R = a*P.X*P.X*P.X
    V = P.Y*P.Y*P.Y
    S = P.Z*P.Z*P.Z
    A = R-V
    A = A*A
    B = R-S
    B = B*B
    C = V-S
    C = C*C
    D = A+C
    E = A+B
    X_new = (R + V + S)*(B + D)
    K = R*C
    L = V*B
    Y_new = d*(K + K - V*(C - E))
    Z_new = d*(L + L - R*(B - D))
    return EC(X_new, Y_new, Z_new)
end

#Repeated tripling
#Input point P, integer e
#Output [3^e]P
function triplePow(P::TwiHesPoint{T}, a::Fp2Elem{T}, d::Fp2Elem{T}, e::Integer) where {T <: Integer}
    Q = P
    for i in 1:e
        Q = triple(Q, a, d)
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

# 3-isogeny evaluation when the kernel is <(0 : -w : 1)>
# Input: Point P, and w such that w^3 = 1, w != 1, parameter a of curve
# Output: Phi(P)
function iso3eval_c1(P::TwiHesPoint{T}, w::Fp2Elem{T}, a::Fp2Elem{T}) where {T <: Integer}
    X, Y, Z = P.X, P.Y, P.Z
    aX3 = a*P.X*P.X*P.X
    Y3 = P.Y*P.Y*P.Y
    Z3 = P.Z*P.Z*P.Z
    X_new = P.X*P.Y*P.Z
    Y_new = aX3 + w*(w*Y3 + Z3)
    Z_new = aX3 + w*(Y3 + w*Z3)
    return EC(X_new, Y_new, Z_new)
end

# 3-isogeny evaluation when the kernel is <(b : -cb : 0)>
# Input: Point P, and cb and b such that c^3 = a
# Output: Phi(P)
function iso3eval_c2(P::TwiHesPoint{T}, cb::Fp2Elem{T}, b::Fp2Elem{T}) where {T <: Integer}
    YZ = P.Y*P.Z
    b2 = b*b
    cb2 = cb*cb
    cbb = cb*b
    X2 = P.X*P.X
    Y2 = P.Y*P.Y
    Z2 = P.Z*P.Z
    b2YZ = b2*YZ
    cb2X2 = cb2*X2
    cbbX = cbb*P.X
    X_new = P.X*b2YZ
    Y_new = cb2X2*P.Z + cbbX*Y2 + b2YZ*P.Z
    Z_new = cb2X2*P.Y + cbbX*Z2 + b2YZ*P.Y
    return EC(X_new, Y_new, Z_new)
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

#Input: A generator point G of order 3^e, and two other points P and Q, curve parameters a and, a non-trivial cuberoot of one, isogeny degree e, and an optimal strategy S
#Output: E/<G>, and phi(P), phi(Q)
function iso3pow_opt(G::TwiHesPoint{T}, P::TwiHesPoint{T}, Q::TwiHesPoint{T}, a::Fp2Elem{T}, d::Fp2Elem{T}, w::Fp2Elem{T}, e::Integer, Strat::Vector{Int64}) where {T <: Integer}
    Pi = P
    Qi = Q
    ai = a
    di = d
    Si = Vector{Tuple{Integer,TwiHesPoint{T}}}()
    push!(Si, (e, G))
    i = 1
    while !isempty(Si)
        (h, Ki) = pop!(Si)
        if h == 1
            b = Ki.X
            if iszero(b)
                #Case 1
                S_new = Vector{Tuple{Integer,TwiHesPoint{T}}}()
                while !isempty(Si)
                    (h, Xi) = popfirst!(Si)
                    push!(S_new, (h-1, iso3eval_c1(Xi, w, ai)))
                end
                Si = S_new
                Pi = iso3eval_c1(Pi, w, ai)
                Qi = iso3eval_c1(Qi, w, ai)
                ai = di*di*di - 27*ai
                di = 3*di
            else
                #Case 2
                cb = -Ki.Y
                if iszero(cb)
                    cb = -Ki.Z
                end
                S_new = Vector{Tuple{Integer,TwiHesPoint{T}}}()
                while !isempty(Si)
                    (h, Xi) = popfirst!(Si)
                    push!(S_new, (h-1, iso3eval_c2(Xi, cb, b)))
                end
                Si = S_new
                Pi = iso3eval_c2(Pi, cb, b)
                Qi = iso3eval_c2(Qi, cb, b)
                c = cb*inv(b)
                dc = di*c
                ai = dc*di + 3*dc*c + 9*ai
                di = di + 6*c
            end
        else
            strat_i = Strat[i]
            if 0 < strat_i < h
                push!(Si, (h, Ki))
                push!(Si, (h-strat_i, triplePow(Ki, ai, di, strat_i)))
                i += 1
            else
                error("Strategy seems to be invalid")
            end
        end
    end
    return Pi, Qi, ai, di
end


# Second part, no other points to carry
#Input: A generator point G of order 3^e, curve parameters a and, a non-trivial cuberoot of one, isogeny degree e, and an optimal strategy S
#Output: E/<G>, and phi(P), phi(Q)
function iso3pow_opt(G::TwiHesPoint{T}, a::Fp2Elem{T}, d::Fp2Elem{T}, w::Fp2Elem{T}, e::Integer, Strat::Vector{Int64}) where {T <: Integer}
    ai = a
    di = d
    Si = Vector{Tuple{Integer,TwiHesPoint{T}}}()
    push!(Si, (e, G))
    i = 1
    while !isempty(Si)
        (h, Ki) = pop!(Si)
        if h == 1
            b = Ki.X
            if iszero(b)
                #Case 1
                S_new = Vector{Tuple{Integer,TwiHesPoint{T}}}()
                while !isempty(Si)
                    (h, Xi) = popfirst!(Si)
                    push!(S_new, (h-1, iso3eval_c1(Xi, w, ai)))
                end
                Si = S_new
                ai = di*di*di - 27*ai
                di = 3*di
            else
                #Case 2
                cb = -Ki.Y
                if iszero(cb)
                    cb = -Ki.Z
                end
                S_new = Vector{Tuple{Integer,TwiHesPoint{T}}}()
                while !isempty(Si)
                    (h, Xi) = popfirst!(Si)
                    push!(S_new, (h-1, iso3eval_c2(Xi, cb, b)))
                end
                Si = S_new
                c = cb*inv(b)
                dc = di*c
                ai = dc*di + 3*dc*c + 9*ai
                di = di + 6*c
            end
        else
            strat_i = Strat[i]
            if 0 < strat_i < h
                push!(Si, (h, Ki))
                push!(Si, (h-strat_i, triplePow(Ki, ai, di, strat_i)))
                i += 1
            else
                error("Strategy seems to be invalid")
            end
        end
    end
    return ai, di
end
