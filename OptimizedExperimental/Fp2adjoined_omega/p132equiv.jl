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
const ea = 63
const eb = 41
const f = big(165)
const p = la^ea*f-1
const F = Fp2(p)

#cuberoot of One
const w = F(big(0),big(1))
if !isone(w^3)
    error("w is not a cuberoot...")
end

#Starting curve: a = 1 for speed, solve for d to get j(E) = some supersingular j-invariant. See sage-script
#Starting j-invariant = 681983795137470565003 + 525699670422529264746*w
#See sage script for construction of d0 from j-Invariant and a0
const a0 = one(F)
const d0 = F(big(51849466119744637087), big(763718193217644493204))

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
PX1 = F(big(177014855531225039999), big(1202619811203149879755))
PY1 = F(big(312237303984514937680), big(1137280334077319004711))
QX1 = F(big(1328653618852880321918), big(671857271963059083134))
QY1 = F(big(905775967463901689844), big(260557785141733845323))
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
