#
# Functions in this file generate simple many-body states
# 

# |↑⟩ = |0⟩ = [1,0], |↓⟩ = |1⟩ = [0,1]   
qubit0::Vector{ComplexF64} = [1, 0]
qubit1::Vector{ComplexF64} = [0, 1]

# Generate a state with maximal bipartite entanglement, and size 2*LA
function maximally_entangled_state(LA::Int64)
    psi = zeros(2^(2*LA))
    max = 2^LA-1
    for i in 0:max
        bit_str = digits(i, base = 2, pad = LA)
        state = 1
        for bit in vcat(bit_str, reverse(bit_str))
            bit == 0 ? state = kron(qubit0, state) : state = kron(qubit1, state)
        end
        psi += reshape(state, (2^(2*LA)))
    end
    psi ./= sqrt(2^LA)
    State(2*LA, psi)
end

function neel_state(L::Int64)
    psi = qubit1
    for r in 2:L
        r % 2 == 0 ? psi = kron(qubit0, psi) : psi = kron(qubit1, psi)
    end
    return State(L, reshape(psi, (2^L)))
end

function anti_neel_state(L::Int64)
    psi = qubit0
    for r in 2:L
        r % 2 == 0 ? psi = kron(qubit1, psi) : psi = kron(qubit0, psi)
    end
    State(L, reshape(psi, (2^L)))
end

function plus_state(L::Int64)
    plus = [1 1]/sqrt(2)
    psi = plus
    for r in 2:L
        psi = kron(plus, psi)
    end
    State(L, reshape(psi, (2^L)))
end

function zero_state(L::Int64)
    psi = qubit0
    for r in 2:L
        psi = kron(qubit0, psi)
    end
    State(L, reshape(psi, (2^L)))
end

# |T⟩
function t_state(L::Int64)
    t = (qubit0 + exp(im*π/4)*qubit1)/sqrt(2)
    psi = t
    for r in 2:L
        psi = kron(t, psi)
    end
    State(L, reshape(psi, (2^L)))
end

function theta_state(L::Int64, theta::Float64)
    local_psi = (qubit0 + exp(im*theta)*qubit1)/sqrt(2)
    psi = local_psi
    for r in 2:L
        psi = kron(local_psi, psi)
    end
    State(L, reshape(psi, (2^L)))
end

function bit_string_state(L::Int64, bit_str::Int64)
    psi = zeros(2^L)
    psi[bit_str] = 1
    return State(L, psi)
end

# |ψ⟩ = |ψ_A⟩ ⊗ |ψ_B⟩
function product_state(stateA::State, stateB::State)
    L = stateA.L + stateB.L
    State(L, reshape(kron(stateB.state, stateA.state), (2^L)))
end

# |Ψ⟩ = |0⟩^(⊗N) ⊗ |T⟩^(⊗P)
function zero_t_state(num_zeros::Int64, num_ts::Int64)
    if(num_zeros == 0)
        return t_state(num_ts)
    elseif(num_ts == 0)
        return zero_state(num_zeros)
    else
        return product_state(zero_state(num_zeros), t_state(num_ts))
    end
end

# |Ψ⟩ = |0⟩^(⊗N0) ⊗ |θ⟩^(⊗Nθ)
function zero_theta_state(num_zeros::Int64, num_thetas::Int64, theta::Float64)
    if(num_zeros == 0)
        return theta_state(num_thetas, theta)
    elseif(num_thetas == 0)
        return zero_state(num_zeros)
    else
        return product_state(zero_state(num_zeros), theta_state(num_thetas, theta))
    end
end

# |Ψ⟩ = |0⟩^(⊗LA) ⊗ |θ⟩^(⊗Nθ) ⊗ |0⟩^(⊗(L - LA - Nθ))
function zero_theta_zero_state(L::Int64, LA::Int64, num_thetas::Int64, theta::Float64)
    if(L - LA - num_thetas == 0)
        return zero_theta_state(LA, num_thetas, theta)
    else
        return product_state(zero_theta_state(LA, num_thetas, theta), zero_state(L - LA - num_thetas))
    end
end

# Random Haar state
function haar_state(L::Int64)
    state = zero_state(L)
    state.state = haar_gate(L) * state.state
    return state
end

function pauli_state(orientation::String)
    basis_str::Vector{String} = ["+x", "-x", "+y", "-y", "+z", "-z"]
    basis = [(qubit0 + qubit1)/sqrt(2), (qubit0 - qubit1)/sqrt(2),
             (qubit0 + im*qubit1)/sqrt(2), (qubit0 - im*qubit1)/sqrt(2),
             qubit0, qubit1]
    index = findfirst(==(orientation), basis_str)
    if(index  == nothing)
        throw(ArgumentError(orientation, "Not a valid spin orientation"))
    end 
    return State(1, basis[index])
end

function pauli_string_state(orientations::String)
    L::Int64 = div(length(orientations), 2)
    psi = pauli_state(orientations[1:2])
    for i in 2:L
        psi = product_state(psi, pauli_state(orientations[2*i-1:2*i]))
    end
    return psi
end
