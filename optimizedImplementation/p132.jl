include("twihes.jl")

############################
#
# Parameters
#
############################

#Prime and field
const la = big(2)
const lb = big(3)
const ea = 63
const eb = 41
const f = big(11)
const p = la^ea*lb^eb*f-1
const F = Fp2(p)

#cuberoot of One
const w = F(big(1850222081870264162797200520152908562431), big(2673953603825912900316077551548884232036))
if !isone(w^3)
    error("w is not a cuberoot...")
end
#Starting curve: a = 1 for speed, solve for d to get j(E) = some supersingular j-invariant. See sage-script
#Starting j-invariant = 2693591125050761171576909111745564865928 + 961716384504775817613723332947597560404*i
#See sage script for construction of d0 from j-Invariant and a0
const a0 = one(F)
const d0 = F(big(422278216911606001622498161204149810012), big(3262366252633517469770463432931360800530))


#Generates torsion basises from random points
function genTorsionBasis(alicePX, alicePY, aliceQX, aliceQY, bobPX, bobPY, bobQX, bobQY)
    Z = one(F)
    aliceP = EC(alicePX, alicePY, Z)
    aliceQ = EC(aliceQX, aliceQY, Z)

    if !isValid(aliceP, a0, d0)
        error("Alice' P is not a valid point!")
    end

    if !isValid(aliceQ, a0, d0)
        error("Alice' Q is not a valid point!")
    end

    Pa = normalized(mul(f*lb^eb, aliceP, a0))
    Qa = normalized(mul(f*lb^eb, aliceQ, a0))

    Patest = doublePow(Pa, a0, ea-1)
    Qatest = doublePow(Qa, a0, ea-1)

    if !isIdentity(double(Patest, a0))
        error("Pa : Something does not look right with the order of the curve...")
    end

    if !isIdentity(double(Qatest, a0))
        error("Qa : Something does not look right with the order of the curve...")
    end

    if isIdentity(Patest)
        error("Alice' point Pa does not have the right order!")
    end
    if isIdentity(Qatest)
        error("Alice' point Qa does not have the right order!")
    end

    if Patest == Qatest
        error("Alice' points are not independent!")
    end

    bobP = EC(bobPX, bobPY, Z)
    bobQ = EC(bobQX, bobQY, Z)

    if !isValid(bobP, a0, d0)
        error("Bob's P is not a valid point!")
    end

    if !isValid(bobQ, a0, d0)
        error("Bob's Q is not a valid point!")
    end

    Pb = normalized(mul(f*la^ea, bobP, a0))
    Qb = normalized(mul(f*la^ea, bobQ, a0))

    Pbtest = triplePow(Pb, a0, d0, eb-1)
    Qbtest = triplePow(Qb, a0, d0, eb-1)

    if !isIdentity(triple(Pbtest, a0, d0))
        error("Pb : Something does not look right with the order of the curve...")
    end

    if !isIdentity(triple(Qbtest, a0, d0))
        error("Qb : Something does not look right with the order of the curve...")
    end

    if isIdentity(Pbtest)
        error("Bob's point Pb does not have the right order!")
    end

    if isIdentity(Qbtest)
        error("Bob's point Qb does not have the right order!")
    end

    if Pbtest == Qbtest
        error("Bob's points are not independent!")
    end

    return Pa, Qa, Pb, Qb
end

#See sage script for generating random points
PX1 = F(big(1805836423159982688664217077118110412720), big(412417461613494901689211800573012059787))
PY1 = F(big(3015485764111398198338164604901199543856), big(2951198480986376391456048710063482208753))
QX1 = F(big(2112425028560317460711754523516298939010), big(3444345531585661127352817874939503284634))
QY1 = F(big(2535178337117861842320927321974655263354), big(1273543743946624092604137125910286830530))
PX2 = F(big(2112666014732165672897417504914718822861), big(2361187632248310541147451912670179879625))
PY2 = F(big(1667563519061147685072119396690356185732), big(2907572659195525987815928592364380130729))
QX2 = F(big(1109723103314231023311003456762614214933), big(1520336319368587385323273075218335709843))
QY2 = F(big(1565635656952058922590166843737964035549), big(3567086747729170057644588319764979489080))
const Pa, Qa, Pb, Qb = genTorsionBasis(PX1, PY1, QX1, QY1, PX2, PY2, QX2, QY2)

# Compute the optimal strategy of size n, given weights p and q
function computeStrat(n, p, q)
    S = Dict{Int64, Array{Int64, 1}}()
    C = Dict{Int64, Int64}()
    S[1] = []
    C[1] = 0
    for i in 2:n
        b = argmin([C[i-b] + C[b] + b*p + (i-b)*q for b in 1:i-1])
        S[i] = append!(append!([b], S[i-b]), S[b])
        C[i] = C[i-b] + C[b] + b*p + (i-b)*q
    end
    return S[n]
end

#The weights are estimated using @benchmark
const strat_A = Vector(computeStrat(ea, 23358, 32132))

const strat_B = Vector(computeStrat(eb, 45620, 34278))
