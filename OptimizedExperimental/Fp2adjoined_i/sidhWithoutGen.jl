include("twihesexp.jl")
using BenchmarkTools

#With strat
function AliceI(s::T, Pa::TwiHesPoint{T}, Qa::TwiHesPoint{T}, a0::Fp2Elem{T}, ea::U, strat_A::Array{Int64}) where {T <: Integer, U <: Integer}
    G = add(Pa, mul(s, Qa, a0), a0)
    a_new, d_new = iso2pow_opt(G, a0, ea, strat_A)
    return a_new, d_new
end

function AliceII(s::T, Pa::TwiHesPoint{T}, Qa::TwiHesPoint{T}, a0::Fp2Elem{T}, ea::Int64, strat_A::Array{Int64}) where {T <: Integer, U <: Integer}
    G = add(Pa, mul(s, Qa, a0), a0)
    a_new, d_new = iso2pow_opt(G, a0, ea, strat_A)
    return jInvariant(a_new, d_new)
end

#Idea: returns the -c parameter of one of S_2, S_3 or S_4 when a=1
function S(i::Int64, w::Fp2Elem{T}) where {T <: Integer}
    return w^(i-2)
end

function addToOracle(OracleList::Array{Fp2Elem{T}}, secretAlice::T, P_A::TwiHesPoint{T}, Q_A::TwiHesPoint{T}, a::Fp2Elem{T}, ea::U, strat_A::Array{Int64}) where {T <: Integer, U <: Integer}
    push!(OracleList, AliceI(secretAlice, P_A, Q_A, a, ea, strat_A)[2])
end

function BobI(list::Array{Int64}, Pa::TwiHesPoint{T}, Qa::TwiHesPoint{T}, a0::Fp2Elem{T}, d0::Fp2Elem{T}, s_a::T, e_a::U, strat_A::Array{Int64}, w::Fp2Elem{T}) where {T <: Integer, U <: Integer}
    Pi = Pa
    Qi = Qa
    di = d0
    ai = a0
    OracleList = Fp2Elem{T}[]
    for i in 1:length(list)
        ci = S(list[i], w)
        Pi = iso3eval_c2(Pi, ci)
        Qi = iso3eval_c2(Qi, ci)
        ai = di*di*ci + 3*di*ci*ci + F(T(9),T(0))
        di = di + 6*ci
        ai = cuberoot(ai)
        Pi = isoHessianEval(Pi, ai)
        Qi = isoHessianEval(Qi, ai)
        di = di*inv(ai)
        addToOracle(OracleList, s_a, Pi, Qi, one(F), e_a, strat_A)
    end
    return Pi, Qi, one(F), di, OracleList
end

function Oracle(dd::Fp2Elem{T}, i::Int64, w::Fp2Elem{T}, OracleList::Array{Fp2Elem{T}}) where {T <: Integer}
    d = OracleList[i]
    for i in 0:2
        if dd == d
            return i
        end
        dd = dd*w
    end
    error("something is not right...")
end

function BobII(list::Array{Int64}, OracleList::Array{Fp2Elem{T}}, a_A::Fp2Elem{T}, d_A::Fp2Elem{T}, w::Fp2Elem{T}) where {T <: Integer}
    di = d_A
    ai = a_A
    for i in 1:length(list)
        ci = S(list[i], w)
        dc = di*ci
        ai = di*dc + 3*dc*ci + F(T(9),T(0))
        di = di + 6*ci
        ai = cuberoot(ai)
        di = di*inv(ai)
        di = di*w^Oracle(di, i, w, OracleList)
    end
    return jInvariant(one(F), di)
end

function transformSBodd(sB::Array{Int64})
    sBnew = Int64[]
    for i in sB
        if i == 3
            push!(sBnew, 4)
        elseif i == 4
            push!(sBnew, 3)
        else
            push!(sBnew, i)
        end
    end
    return sBnew
end

function example()
    sA = rand(big(1):big(p))
    sB = [rand(2:4) for i in 1:eb]
    println("\nAlice' secret: $sA")
    println("\nBob's secret: $sB")
    a_A, d_A = AliceI(sA, Pa, Qa, a0, ea, strat_A)
    P_A, Q_A, a_B, d_B, OracleList = BobI(sB, Pa, Qa, a0, d0, sA, ea, strat_A, w)

    println("\nAlice sends:")
    println("E_A : ($(a_A))X^3 + Y^3 + Z^3 = ($(d_A))XYZ\n")
    println("Bob sends:")
    println("P_A = $P_A")
    println("Q_A = $Q_A")
    println("E_B : ($(a_B))X^3 + Y^3 + Z^3 = ($(d_B))XYZ\n")
    println("--------------------------------\n")

    if ea % 2 == 1
        sB = transformSBodd(sB)
    end

    key_a = AliceII(sA, P_A, Q_A, a0, ea, strat_A)
    key_b = BobII(sB, OracleList, a_A, d_A, w)

    println("Results:")
    println("key_a : $key_a")
    println("key_b : $key_b")
end

function benchmarkBob()
    sA = rand(big(1):big(p))
    sB = [rand(2:4) for i in 1:eb]
    a_A, d_A = AliceI(sA, Pa, Qa, a0, ea, strat_A)
    P_A, Q_A, a_B, d_B, OracleList = BobI(sB, Pa, Qa, a0, d0, sA, ea, strat_A, w)

    if ea % 2 == 1
        sB = transformSBodd(sB)
    end

    key_a = AliceII(sA, P_A, Q_A, a0, ea, strat_A)
    key_b = BobII(sB, OracleList, a_A, d_A, w)

    if key_a != key_b
        error("Shared secret was not the same")
    end

    return @benchmark BobII($sB, $OracleList, $a_A, $d_A, $w)
end

function benchmarkAlice()
    sA = rand(big(1):big(p))
    sB = [rand(2:4) for i in 1:eb]
    a_A, d_A = AliceI(sA, Pa, Qa, a0, ea, strat_A)
    P_A, Q_A, a_B, d_B, OracleList = BobI(sB, Pa, Qa, a0, d0, sA, ea, strat_A, w)

    if ea % 2 == 1
        sB = transformSBodd(sB)
    end

    key_a = AliceII(sA, P_A, Q_A, a0, ea, strat_A)
    key_b = BobII(sB, OracleList, a_A, d_A, w)

    if key_a != key_b
        error("Shared secret was not the same")
    end

    return @benchmark AliceII($sA, $P_A, $Q_A, $a0, $ea, $strat_A)
end


function main()
    choice = "0"
    while !(choice in ["1","2","3","4","5"])
        println("Which parameter-sets do you want?")
        println("[1] p132, [2] p434, [3] p503, [4] p610, [5] p751")
        choice = readline()
        if choice == "1"
            include("p132equiv.jl")
        elseif choice == "2"
            include("p434equiv.jl")
        elseif choice == "3"
            include("p503equiv.jl")
        elseif choice == "4"
            include("p610equiv.jl")
        elseif choice == "5"
            include("p751equiv.jl")
        else
            println("pick a number between 1 and 5")
            println("----------------------------------------------------------------")
        end
    end
    choice = "0"
    while choice != "4"
        println("What do you want to do?")
        println("[1] Run as example, [2] Benchmark Alice, [3] Benchmark Bob, [4] Close")
        choice = readline()
        if choice == "1"
            example()
        elseif choice == "2"
            println("\n")
            show(stdout,MIME"text/plain"(),benchmarkAlice())
            println("\n")
        elseif choice == "3"
            println("\n")
            show(stdout,MIME"text/plain"(),benchmarkBob())
            println("\n")
        elseif choice == "4"
            return
        else
            println("pick a number between 1 and 4")
            println("----------------------------------------------------------------")
        end
    end
end

main()
