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
const w = F(big(27084430341343327003025253873091092029340323762182288388458123414615853761653324845414451511295),big(43469381832290841702031770221732784416315090589113876447022926993991394359825475508483569005976))
if !isone(w^3)
    error("w is not a cuberoot...")
end

#Starting curve: a = 1 for speed, solve for d to get j(E) = some supersingular j-invariant. See sage-script
#Starting j-invariant = 33893307480695892084930852953913044770858278533391760625844466115645790391283264792416646998918 + 49893416820742440067790454993290701172148055052125782460670673388722101266446998614247240705282*i
#See sage script for construction of d0 from j-Invariant and a0
const a0 = one(F)
const d0 = F(big(9582897442711434728403124384511139620628929294208952482427567967616014771315321094720418930686), big(31054706339263183158896331073932276558215085803835343090123395621182233568259675682838496316922))

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
PX1 = F(big(24987478090064704368145591133367866412348882360494478779537919087897415026700920919820345306964), big(45156443088543560408882750447981258295501394838557271722956534941967060115719852316788629482718))
PY1 = F(big(2170702476518712889664320898445979081591013738410417670804682017474321794329043063458007488306), big(49082543204777057657712331797698773742733613316359103302811716027579956455981335192338668013907))
QX1 = F(big(53528141194563294088405111131047362264395244595521304890182543770751518491009651832303668435420), big(13652615373534279301709419661007866284982621832443381365418404541743107498645664796623590178471))
QY1 = F(big(18909889937062329608189290706836126453204719304330955995617636788633998498984933840787167162794), big(16073651271491968770496512250315521100829891818876969764255765580174478538389443332128361994708))
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
