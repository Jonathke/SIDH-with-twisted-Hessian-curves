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
const ea = 305
const eb = 192
const f = big(831)
const p = la^ea*f-1
const F = Fp2(p)

#cuberoot of One
const w = F(big(0),big(1))
if !isone(w^3)
    error("w is not a cuberoot...")
end

#Starting curve: a = 1 for speed, solve for d to get j(E) = some supersingular j-invariant. See sage-script
#Starting j-invariant =6524128774462177850216699667974953090043279650578252070684471300279639234307848380131542810376 + 39465243035713526479017432822883889983669536408709373552072994697416033615854320664321601455920*w
#See sage script for construction of d0 from j-Invariant and a0
const a0 = one(F)
const d0 = F(big(1431317826531827479268470564041495987770387357124362784334083726662800734151339932879547105456), big(51664592964949788114266166910073594440142572679665527046679550053192549887213743832589932028458))

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
PX1 = F(big(6392336050481667159825927734381233714892801816974840185362818522638755097934850968664562994071), big(9964681269262279261764151625222066351738951262982552259454137779293345367974100077397912408993))
PY1 = F(big(34352663371017071108764651376057085936969514304846959295476043951330873473886162051626461572802), big(13948232429937690716941033492581256276040107866654060300309241461719735705535322736545024524583))
QX1 = F(big(3232339357389555004070339343093561631492892383570542817704076142493998627625650457561961732458), big(36415596284512728235124093054409146315471110438123034852791646213545858129087876953033466918641))
QY1 = F(big(50802091744828683992108894196984444634061617992673604106895466802127033503494503543639753289083), big(41086673058833284235667898049886726817522123989298391499252687503009772952336376360337878835002))
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
