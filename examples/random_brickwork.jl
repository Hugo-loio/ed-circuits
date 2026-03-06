import EDCircuit as ed

#=
This example serves as a benchmark of the performance of this package.
It compute the runtime of a brickwork evolution of depth D 
of random Haar 2-qubit gates on a system of N qubits 
(each layer acts both on odd and even neighboring sites) 
=#

N = 18 # Total number of qubits in a 1D chain
D = N # Circuit depth
pbc = true # Periodic boundary conditions
 
# First order trotterization scheme with 2-qubit gates
circuit = ed.Circuit(pbc, N, 0)
trotter_sites = ed.trotter_sites_vector(circuit, 2)  

state = ed.zero_state(N) # Initial state |0>^{⊗N}

println("\nStarting benchmark for N = $N, D = $D...")

# Timing the evolution loop
runtime = @elapsed begin
    for l in 1:D
        ed.apply_random_haar_gates!(state, trotter_sites)
    end
end

println("Total Runtime: $(round(runtime, digits=4)) seconds")
println("Average Time per Layer: $(round(runtime/D, digits=6)) seconds\n")
