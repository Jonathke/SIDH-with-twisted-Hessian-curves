include("TwistedHessianEC.jl")
using .TwistedHessianEC

function partI(m::W, n:: W, P::S, Q::S, Pother::S, Qother::S, l::U, e::V, E::T) where {T <: EllipticCurve, S <: EllipticCurveElem, U <: Integer, V <: Integer, W <: Integer}
    K = P*m + Q*n
    Ei = E
    Pi = Pother
    Qi = Qother
    for i in V(1):e
        Ki = (l^(e-i))*K
        phi = Isogeny(Ei, Ki)
        K = phi(K)
        Pi = phi(Pi)
        Qi = phi(Qi)
        Ei = phi.codomain
    end
    return Ei, Pi, Qi
end

function partII(m::W, n:: W, P::S, Q::S, l::U, e::V, E::T) where {T <: EllipticCurve, S <: EllipticCurveElem, U <: Integer, V <: Integer, W <: Integer}
    K = P*m + Q*n
    Ei = E
    for i in V(1):e
        Ki = (l^(e-i))*K
        phi = Isogeny(Ei, Ki)
        K = phi(K)
        Ei = phi.codomain
    end
    return jInvariant(Ei)
end

# Example:
#Public params
la = big(2)
lb = big(3)
ea = 63
eb = 41
OGp = la^ea*lb^eb*11-1
F = Fpsquare(OGp)
E = TwiHesEC(F, one(F), F([3, 2458499195747164226075539892235780232099]))

XP = F([big(264972986297324027107883370923360079519), big(483679027370140600263719808627289108176)])
YP = F([big(2550992266319579452440794756588926041584), big(2214394007774501716861032109362868970653)])
Z = one(F)
P = E(XP, YP, Z)

XQ = F([big(1287714346072412371952478415842888698876), big(427368407186067938537508979187813002488)])
YQ = F([big(528075014001520129686916245492679497447), big(857779714716146821676003610752395199871)])
Q = E(XQ, YQ, Z)

Pa = normalized(P*(big(3)^41*11))
Qa = normalized(Q*(big(3)^41*11))

Pb = normalized(P*(big(2)^63*11))
Qb = normalized(Q*(big(2)^63*11))


# Alice and Bob chooses secret params
# They cannot both be divisible prime by "their" prime...
function genSecret(la, lb)
    ma = rand(big(1):OGp-1)
    na = rand(big(1):OGp-1)
    while ma%la == 0 && na%la == 0
        ma = rand(big(1):OGp-1)
        na = rand(big(1):OGp-1)
    end
    mb = rand(big(1):OGp-1)
    nb = rand(big(1):OGp-1)
    while mb%lb == 0 && nb%lb == 0
        mb = rand(big(1):OGp-1)
        nb = rand(big(1):OGp-1)
    end
    return ma, na, mb, nb
end

ma, na, mb, nb = genSecret(la, lb)


# Alice do her part
E_A, P_B, Q_B = partI(ma, na, Pa, Qa, Pb, Qb, la, ea, E)
# Bob does his part
E_B, P_A, Q_A = partI(mb, nb, Pb, Qb, Pa, Qa, lb, eb, E)

#Both derive the shared key, hopefully it is equal
key_a = partII(ma, na, P_A, Q_A, la, ea, E_B)
key_b = partII(mb, nb, P_B, Q_B, lb, eb, E_A)
println("Alice key: $key_a")
println("Bob's key: $key_b")


# Del
