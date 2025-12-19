# 
# This file defines some useful gates
# 

# Random haar n-qubit gate
function haar_gate(num_sites::Int64)
    dim = 2^num_sites
    z = (randn(dim, dim) + im*randn(dim, dim)) / sqrt(2.0)
    q, r = qr(z)
    d = Diagonal([sign(x) for x in diag(r)])
    return Matrix(q) * d
end

# Random haar 2-qubit gate
function haar_gate_2qb()
    return haar_gate(2)
end


# Random complex phase
function random_complex_phase()
    phase = rand(Complex{Float64})
    return phase/norm(phase)
end

# 2-qubit haar U1 preserving gate
function haar_U1_gate_2qb()
    z = (randn(2,2) + im*randn(2,2)) / sqrt(2.0)
    q, r = qr(z)
    d = Diagonal([sign(x) for x in diag(r)])
    gate = zeros(ComplexF64, 4,4)
    gate[1,1] = random_complex_phase()
    gate[4,4] = random_complex_phase()
    gate[2:3,2:3] = Matrix(q) * d
    return gate
end

# T gate 
tgate = zeros(ComplexF64, 2, 2)
tgate[1,1] = 1
tgate[2,2] = exp(im*π/4)

function phase_gate(θ::Float64)
    return [1 0 ; 0 exp(im*θ)]
end
