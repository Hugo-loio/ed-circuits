using InteractiveUtils
using BenchmarkTools
using LinearAlgebra
using StatsBase
import EDCircuit as ed
 
L = 10
state = ed.neel_state(L)

#@code_warntype ed.haar_gate(2)
#@code_warntype ed.apply_gate!(state, [4,6], ed.haar_gate(2))
#@code_warntype ed.measure!(state, 10)
#@code_warntype ed.exact_projected_ensemble(state, [1,2])

circuit = ed.Circuit(false, L, 0)
trotter_sites = ed.trotter_sites_vector(circuit, 2)
for _ in 1:L
    ed.apply_random_haar_gates!(state, trotter_sites)
end

@code_warntype ed.rand_proj_group_overlaps_frame_potential(state, [1], 6, [10,100])

@code_warntype ed.exact_proj_group_overlaps_frame_potential(state, [1], 6, [10,100])

