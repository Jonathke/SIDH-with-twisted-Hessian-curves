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
const w = F(big(0),big(1))
if !isone(w^3)
    error("w is not a cuberoot...")
end

#Starting curve: a = 1 for speed, solve for d to get j(E) = some supersingular j-invariant. See sage-script
#Starting j-invariant = 386702733408549460782724268732857992185964070161051245286529361669022244913896 + 580419299742272852086148469132779561091771935013448986984506353121177048771958*w
#See sage script for construction of d0 from j-Invariant and a0
const a0 = one(F)
const d0 = F(big(17295511272188241945560521111123569448987393250708156671063629469542861050669), big(349781430796024944336522884354347512620320027511390987651460415007763363894936))

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
PX1 = F(big(617062450513364498991242959109928480156747435503827845620272706880260818054087), big(210997965445286911603227942566208008681285659915304257006978481066551090223362))
PY1 = F(big(492371184480695605318885506056597828272549417762827810563525060531233108248103), big(92900268813052460312972905161689196944082919281989013977250198785826086779770))
QX1 = F(big(127833779941631013946288599048105980912964446225551112260445274440126767802139), big(102229667146288756900372263903576782020907595160379577145862433826866036173834))
QY1 = F(big(307503179749308467456513487591472521505681708512696786529178002999800233174999), big(52634220846657808547056742322019734562168630750820301350085458668631624844436))
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
