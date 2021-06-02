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
const ea = 372
const eb = 239
const f = big(933)
const p = la^ea*f-1
const F = Fp2(p)

#cuberoot of One
const w = F(big(0),big(1))
if !isone(w^3)
    error("w is not a cuberoot...")
end

#Starting curve: a = 1 for speed, solve for d to get j(E) = some supersingular j-invariant. See sage-script
#Starting j-invariant = 8847223067209440401448229588737313272559879717379846897699863504053610128797852792001224868454241111711674404439151 + 6504793230146749584184744164990186194831639488939465878632605461025444666200467465185374670369980430549746221725605*w
#See sage script for construction of d0 from j-Invariant and a0
const a0 = one(F)
const d0 = F(big(73556440198263896216426411979837259119151083503386436972423522873961150562697915691134122217025043180758676456962), big(7365407801211656443489327267965348700482332837331946019398173344961922773597286373170842211746143557209138074396819))

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
PX1 = F(big(7552506659735093889758248207443140101485815964709456410994626916033548008217670936854954221019153609689168549709923), big(3127295300597838738635529884165612772164775098253635177650537509490946653278125103277078936770129940403886285748008))
PY1 = F(big(4335080515298023914789642233386072822783819110885463188275452087633232775973249616247740473168837899241267585477136), big(8511056511339186759429385084799128950140512962083879087521472255203235735603955755091104392619026040526119505515430))
QX1 = F(big(885439283334681463147952303783777282418326871195947052828393343356395115437437294785388621396640008753002437611419), big(683601484854071613822828240177114891640408423460914316412803470427687975251890723100714427278969392553887593851202))
QY1 = F(big(2159582834891235194701010380072949557378667671420445189111379944801618363317225068688477014536560727200433992201094), big(4009897625395987367565414085892970785531526625836881837151754666888765024041657855946948803354989404799753309872286))
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
