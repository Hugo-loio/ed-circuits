import EDCircuit as ed

#=
In this example we consider the many-body dynamics of a brickwork circuit of 
random 2-qubit gates with mid-circuit measurments with probability p. 
Starting from a simple product state you can see how the bipartite entanglement 
entropy grows linearly in time until saturation at either area-law or
volume-law, depending on the value of p.
=#

function print_bip_S(state::ed.State, layer::Int64)
    println("Bipartite entanglement entropy at layer ", layer,
            ": ", ed.bip_ent_entropy(state, sites_bip))
end

L = 18 # Total number of qubits in a 1D chain
Lbip = div(L,2) # Number of qubits in a bipartition
pbc = true # Periodic boundary conditions
D = L # Circuit depth
p = 0.1 # Measurement probability

sites_bip = collect(1:Lbip)

# First order trotterization scheme with 2-qubit gates
circuit = ed.Circuit(pbc, L, 0)
trotter_sites = ed.trotter_sites_vector(circuit, 2)  


for p in [0.1, 0.6]
    # Initial state |0>^{⊗L}
    state = ed.zero_state(L)

    println("\nMeasurement probability ", p)
    print_bip_S(state, 0)
    for l in 1:D
        ed.apply_random_haar_gates!(state, trotter_sites)
        ed.measure_layer!(circuit, state, p)
        print_bip_S(state, l)
    end
end

