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
const w = F(big(4487557590482916150519573292655516562517016203533235075829560271748200489844591980035340626929641881454172841181183),big(5738677033203073227609718885528424252917843825926148651730561078873371523986093711547250640357294362739990427211904))
if !isone(w^3)
    error("w is not a cuberoot...")
end

#Starting curve: a = 1 for speed, solve for d to get j(E) = some supersingular j-invariant. See sage-script
#Starting j-invariant = 5078601030025936612951094422122457493757846004081090863385038905475160431769027333596329668623319709636198670796335 + 730435898302185493450854427674126914241218266865584295975053215582853664936059337902183419108185317207844682910681*i
#See sage script for construction of d0 from j-Invariant and a0
const a0 = one(F)
const d0 = F(big(215204150778632632067183850681783022848182107473406273571526638557875768196798935580442884573473650069368123798743), big(4395200443804932030300697937166824270064850253643806210095266826983082100152405861232785273741882382287841446299456))

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
PX1 = F(big(881950242471972311025635055360385496242977735496495368465907993398961143296737488051828563879732739480367323983878), big(2760196959225527514449993758945154999464866359809639723333677922497968337828709066990814940672436711491108697193597))
PY1 = F(big(2866394269112510555560585835799329040580044781434580664844320118821746139039440025584732107770883210318523356600665), big(7745137799217928158068986469504819438862624704323536011280272734488091546758794307191543146806559301552218995964710))
QX1 = F(big(6099686467544935488250746559066884528337653299859103018409795832698827610498644172702196448512614568599661131423340), big(5071747096883336249394454607847123989533962100221254541289002291417576182986764996156874663244863045274511842996263))
QY1 = F(big(2048487070989477784952038398019216517831515602028577702508207526150192363199729101585599614700610074488676246724945), big(497973770723935457779837924788248583895180744595801874111571202757285894627776214594333839715407911872990454725884))
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
