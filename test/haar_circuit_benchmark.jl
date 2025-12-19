using InteractiveUtils
using BenchmarkTools
using LinearAlgebra
using StatsBase
import EDCircuit as ed
using Profile
using FileIO

L::Int64 = 24
depth::Int64 = L

function run_circuit(depth, verbose=false)
    state = ed.neel_state(L)
    circuit = ed.Circuit(true, L, 0)
    sites = ed.trotter_sites_vector(circuit, 2)
    for step in 1:depth
        if(verbose)
            println(step)
        end
        ed.apply_random_haar_gates!(state, sites)
    end
end
run_circuit(1)

#Profile.clear()
#@profile run_circuit(depth)
#save("haar_circuit.jlprof",  Profile.retrieve()...)

@time run_circuit(depth, true)
