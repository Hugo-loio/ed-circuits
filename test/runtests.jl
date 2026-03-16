#using Pkg
#Pkg.add("Combinatorics")

using Combinatorics
using Test

import EDCircuit as ed


function pauli_string(string::String)
    sigmas = [[0 1; 1 0], [0 -1im; 1im 0], [1 0; 0 -1]]
    res = 1
    for pauli in string
        pauli = lowercase(pauli)
        if pauli == 'x'
            res = kron(res, sigmas[1])
        elseif pauli == 'y'
            res = kron(res, sigmas[2])
        elseif pauli == 'z'
            res = kron(res, sigmas[3])
        else
            error("Pauli operator '$pauli' not recognized!")
        end 
    end
    return Matrix{ComplexF64}(res)
end

function test_string_ev(base::ed.State, string::String, sites::Vector{Int64})
    state1 = deepcopy(base)
    ed.apply_gate!(state1, sites, pauli_string(string))
    eval1 = ed.overlap(base, state1)
    eval2 = ed.pauli_string_ev(base, string, sites)
    return abs(eval1 - eval2) < 1E-10
end


@testset "EDCircuit.jl" begin
    L = 5
    paulis = ['x', 'y', 'z']
    base = ed.haar_state(L)

    for len in 1:L
        strings = [join(x) for x in Iterators.product(fill(paulis, len)...)]
        sitess = collect(permutations(collect(1:L), len))
        for sites in sitess, string in strings
            @test test_string_ev(base, string, sites)
        end
    end
end
