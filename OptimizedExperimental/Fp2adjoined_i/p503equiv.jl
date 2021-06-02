include("twihesexp.jl")

############################
#
# Parameters
#
############################

#Searches for a suitable f
function findf(la, ea, bound)
    for i in 1:bound
        p = la^ea*i*3-1
        if p%9 == 5 && isprime(p)
            return i*3
        end
    end
end

#Prime and field
const la = big(2)
const ea = 250
const eb = 159
const f = big(366)
const p = la^ea*f-1
const F = Fp2(p)

#cuberoot of One
const w = F(big(331093005162950996289273285259216986517943862403315987800324029272626605064191),big(104623851773843360140102269683274291590290222580784450012157868666242955860563))
if !isone(w^3)
    error("w is not a cuberoot...")
end

#Starting curve: a = 1 for speed, solve for d to get j(E) = some supersingular j-invariant. See sage-script
#Starting j-invariant = 368293176498927158256933784097060604462499111768855219774436605057193395479918 + 525122908409866890140632891334529838028770175596997650212099483929576712609501*i
#See sage script for construction of d0 from j-Invariant and a0
const a0 = one(F)
const d0 = F(big(87443932825301580207118065997125411432163462135340894657316799133032834336231), big(502682738134202325970433028490156695512656145794126134896470980539619709732022))

#Generates torsion basises from random points
function genTorsionBasis(alicePX, alicePY, aliceQX, aliceQY)
    Z = one(F)
    aliceP = EC(alicePX, alicePY, Z)
    aliceQ = EC(aliceQX, aliceQY, Z)

    if !isValid(aliceP, a0, d0)
        error("Alice' P is not a valid point!")
    end

    if !isValid(aliceQ, a0, d0)
        error("Alice' Q is not a valid point!")
    end

    Pa = normalized(mul(f, aliceP, a0))
    Qa = normalized(mul(f, aliceQ, a0))

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

    return Pa, Qa
end

#See sage script for generating random points
PX1 = F(big(613744005296302797246819685030604636246204864533350040851263459468697566136457), big(501561563001396309984769357114178587853450197878880784473010702997914231623319))
PY1 = F(big(60795157551561278712985108503724765174093818495989924665590590907821333813183), big(483603975005329087586603834861406666980912662071633613729547534589142805135870))
QX1 = F(big(200516930093991042592609952714723907361941616540955914517379474424255350232393), big(436553646898939476886856239929834413093735295319263029998454094826078440937187))
QY1 = F(big(488228688869854242488382183909375596955556471268393454554170488930576353755991), big(388865888441358174971910959705192289058933214286718079800261178112284469197141))
const Pa, Qa = genTorsionBasis(PX1, PY1, QX1, QY1)

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
