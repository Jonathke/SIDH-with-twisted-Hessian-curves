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
const ea = 216
const eb = 137
const f = big(879)
const p = la^ea*f-1
const F = Fp2(p)

#cuberoot of One
const w = F(big(46284752188330883553734973166973155004154394428541334396269674627071),big(23404647333901950368279213242774107664334304380270304764970672375197))
if !isone(w^3)
    error("w is not a cuberoot...")
end

#Starting curve: a = 1 for speed, solve for d to get j(E) = some supersingular j-invariant. See sage-script
#Starting j-invariant = 14064682111529907679862563344501303303796260730204345953395685638569 + 76096427302910228682162677922758597934852762832035895265033628093994*i
#See sage script for construction of d0 from j-Invariant and a0
const a0 = one(F)
const d0 = F(big(43724639639236444273424904448372818046331617976508621957630798715), big(29708227635556942990944295980410771388696051411388522968285005956166))

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
PX1 = F(big(16633178152203574572039013028007631092696491990507931552466794688583), big(37012916792942429966222346695142422448749627678575350428167984878958))
PY1 = F(big(57022629875491273926838985745691278822530512869563487052868582083704), big(19933533038423915493571572118202112239770435037617206937438373600403))
QX1 = F(big(43618056518659814000877608901074696861003674853424119434042477048020), big(39158734920895635010661994810826191162027597218353879351858614361453))
QY1 = F(big(60772941628506071727569214606909953796153806879412032951445411916955), big(16149854102456248625790194427351319162676700476190050641406454011899))
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
