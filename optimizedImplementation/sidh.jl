include("twihes.jl")
#include("pxxx.jl") <- Globaly defines constant parameters. Useful for when running in REPL
using BenchmarkTools

#Alice' first part of the key-exchange
function ephimeral_A(s::T, Pa::TwiHesPoint{T}, Qa::TwiHesPoint{T}, Pb::TwiHesPoint{T}, Qb::TwiHesPoint{T}, a0::Fp2Elem{T}, ea::U, strat_A::Array{Int64}) where {T <: Integer, U <:Integer}
    G = add(Pa, mul(s, Qa, a0), a0)
    P_B, Q_B, a_new = iso2pow_opt(G, Pb, Qb, a0, ea, strat_A)
    d_new = recover_d(P_B, a_new) #Ignores d-parameter during calculation, so must recover it at the end
    return P_B, Q_B, a_new, d_new
end

#Bob's first part of the key-exchange
function ephimeral_B(s::T, Pb::TwiHesPoint{T}, Qb::TwiHesPoint{T}, Pa::TwiHesPoint{T}, Qa::TwiHesPoint{T}, a0::Fp2Elem{T}, d0::Fp2Elem{T}, eb::U, strat_B::Array{Int64}, w::Fp2Elem{T}) where {T <: Integer, U <:Integer}
    G = add(Pb, mul(s, Qb, a0), a0)
    P_A, Q_A, a_new, d_new = iso3pow_opt(G, Pa, Qa, a0, d0, w, eb, strat_B)
    return P_A, Q_A, a_new, d_new
end

#Alice' second part of the key-exchange
function secretAgreement_A(s::Integer, P_A::TwiHesPoint{T}, Q_A::TwiHesPoint{T}, a_B::Fp2Elem{T}, ea::U, strat_A::Array{Int64}) where {T <: Integer, U <:Integer}
    G = add(P_A, mul(s, Q_A, a_B), a_B)
    a_new, d_new = iso2pow_opt(G, a_B, ea, strat_A) #Recovering d is build in the function, see twihes.jl
    return jInvariant(a_new, d_new)
end

#Bob's second part of the key-exchange
function secretAgreement_B(s::Integer, P_B::TwiHesPoint{T}, Q_B::TwiHesPoint{T}, a_A::Fp2Elem{T}, d_A::Fp2Elem{T}, eb::U, strat_B::Array{Int64}, w::Fp2Elem{T}) where {T <: Integer, U <:Integer}
    G = add(P_B, mul(s, Q_B, a_A), a_A)
    a_new, d_new = iso3pow_opt(G, a_A, d_A, w, eb, strat_B)
    return jInvariant(a_new, d_new)
end

#Full key-exhange example
function example()
    sA = rand(1:la^ea)
    sB = rand(1:lb^eb)
    println("\nAlice' secret: $sA")
    println("\nBob's secret: $sB")

    P_B, Q_B, a_A, d_A = ephimeral_A(sA, Pa, Qa, Pb, Qb, a0, ea, strat_A)
    P_A, Q_A, a_B, d_B = ephimeral_B(sB, Pb, Qb, Pa, Qa, a0, d0, eb, strat_B, w)

    println("\nAlice sends:")
    println("phi_A(P_B) = $P_B")
    println("phi_A(Q_B) = $Q_B")
    println("E_A : ($(a_A))X^3 + Y^3 + Z^3 = ($(d_A))XYZ\n")
    println("Bob sends:")
    println("phi_B(P_A) = $P_A")
    println("phi_B(Q_A) = $Q_A")
    println("E_B : ($(a_B))X^3 + Y^3 + Z^3 = ($(d_B))XYZ\n")
    println("--------------------------------\n")

    j_A = secretAgreement_A(sA, P_A, Q_A, a_B, ea, strat_A)
    j_B = secretAgreement_B(sB, P_B, Q_B, a_A, d_A, eb, strat_B, w)

    println("Results:")
    println("key_A: $j_A")
    println("key_B: $j_B \n \n")
end

#Timing the second part of Bob's key-exhange to compare with OptimizedExperimental
function benchmarkBob()
    sA = rand(1:la^ea)
    sB = rand(1:lb^eb)
    P_B, Q_B, a_A, d_A = ephimeral_A(sA, Pa, Qa, Pb, Qb, a0, ea, strat_A)
    P_A, Q_A, a_B, d_B = ephimeral_B(sB, Pb, Qb, Pa, Qa, a0, d0, eb, strat_B, w)
    j_A = secretAgreement_A(sA, P_A, Q_A, a_B, ea, strat_A)
    j_B = secretAgreement_B(sB, P_B, Q_B, a_A, d_A, eb, strat_B, w)
    if j_A != j_B
        error("Shared secret was not the same")
    end

    return @benchmark secretAgreement_B($sB, $P_B, $Q_B, $a_A, $d_A, $eb, $strat_B, $w)
end

#Timing the second part of Alice' key-exhange to compare with OptimizedExperimental
function benchmarkAlice()
    sA = rand(1:la^ea)
    sB = rand(1:lb^eb)
    P_B, Q_B, a_A, d_A = ephimeral_A(sA, Pa, Qa, Pb, Qb, a0, ea, strat_A)
    P_A, Q_A, a_B, d_B = ephimeral_B(sB, Pb, Qb, Pa, Qa, a0, d0, eb, strat_B, w)
    j_A = secretAgreement_A(sA, P_A, Q_A, a_B, ea, strat_A)
    j_B = secretAgreement_B(sB, P_B, Q_B, a_A, d_A, eb, strat_B, w)
    if j_A != j_B
        error("Shared secret was not the same")
    end

    return @benchmark secretAgreement_A($sA, $P_A, $Q_A, $a_B, $eb, $strat_B)
end

function main()
    choice = "0"
    while !(choice in ["1","2","3","4","5"])
        println("Which parameter-sets do you want?")
        println("[1] p132, [2] p434, [3] p503, [4] p610, [5] p751")
        choice = readline()
        if choice == "1"
            include("p132.jl")
        elseif choice == "2"
            include("p434.jl")
        elseif choice == "3"
            include("p503.jl")
        elseif choice == "4"
            include("p610.jl")
        elseif choice == "5"
            include("p751.jl")
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
