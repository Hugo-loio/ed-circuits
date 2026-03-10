using BenchmarkTools

import EDCircuit as ed

#=
This example serves as a benchmark of the performance of this package.
It simulates the brickwork evolution of depth D 
of random Haar 2-qubit gates on a system of N qubits 
(each layer acts both on odd and even neighboring sites).
It also benchmarks hybrid evolutation with an additional measurement layer.
=#

N = 18 # Total number of qubits in a 1D chain
D = N # Circuit depth
pbc = true # Periodic boundary conditions
p = 0.5 # Measurement probability
 
# First order trotterization scheme with 2-qubit gates
circuit = ed.Circuit(pbc, N, 0)
trotter_sites = ed.trotter_sites_vector(circuit, 2)  

state = ed.zero_state(N) # Initial state |0>^{⊗N}

println("\nStarting unitary benchmark for N = $N, D = $D...")
benchmark = @benchmark begin
    for l in 1:D
        ed.apply_random_haar_gates!(state, trotter_sites)
    end
end
display(benchmark)

println("\nStarting hybrid benchmark for N = $N, D = $D, p = $p...")
benchmark = @benchmark begin
    for l in 1:D
        ed.apply_random_haar_gates!(state, trotter_sites)
        ed.measure_layer!(circuit, state, p)
    end
end
display(benchmark)
