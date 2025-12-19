import EDCircuit as ed

#=
In this example we consider the many-body dynamics of a brickwork circuit of 
random 2-qubit gates. Starting from a simple product state you can see how 
the bipartite entanglement entropy grows linearly in time until saturation.
Furthermore, you can also see how local ensembles converge to the Haar
ensemble exponentially fast, i.e. the system deep thermalizes.
=#

function print_bip_S(state::ed.State)
    println("Bipartite entanglement entropy ",
            ed.bip_ent_entropy(state, sites_bip))
end

function print_haar_dist(state::ed.State)
    frames = ed.exact_proj_overlaps_frame_potential(state, sitesA, kmax)
    dists = frames./haar_frames .- 1
    println("Local distance to the Haar ensemble")
    for k in 1:kmax 
        println("\tk = ", k, ": ", dists[k])
    end
end

L = 16 # Total number of qubits in a 1D chain
Lbip = div(L,2) # Number of qubits in a bipartition
LA = 2 # Number of qubits in a left partition A
pbc = true # Periodic boundary conditions
D = L # Circuit depth
kmax = 4 # Max order Haar distance

sites_bip = collect(1:Lbip)
sitesA = collect(1:LA)
dA = 2^LA # Hilbert space dimension of A

haar_frames = [ed.haar_frame_potential(dA, k) for k in 1:kmax]

# First order trotterization scheme with 2-qubit gates
circuit = ed.Circuit(pbc, L, 0)
trotter_sites = ed.trotter_sites_vector(circuit, 2)  

# Initial state |0>^{⊗L}
state = ed.zero_state(L)

println("Initial state")
print_bip_S(state)
print_haar_dist(state)
for l in 1:D
    ed.apply_random_haar_gates!(state, trotter_sites)
    println("\nLayer ", l, ":")
    print_bip_S(state)
    print_haar_dist(state)
end

