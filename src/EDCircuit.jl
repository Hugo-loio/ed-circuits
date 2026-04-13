#
# This module implements exact diagonalization simulations for 1D spin-1/2 chains
#

module EDCircuit

using Dates
using LinearAlgebra
using SparseArrays
using SpecialFunctions: factorial
using Combinatorics: permutations

# Container for circuit properties
# Often the struct is not needed - TODO: remove in the future
mutable struct Circuit
    pbc::Bool
    L::Int64
    L_ancillae::Int64
    oddfirst::Bool
    sites::NTuple
    oddsites::NTuple
    evensites::NTuple
    function Circuit(pbc, L, L_ancillae; oddfirst=true)
        sites = tuple(collect((1+L_ancillae):(L+L_ancillae))...)
        oddsites = tuple(collect((1+L_ancillae):2:(L+L_ancillae))...)
        evensites = tuple(collect((2+L_ancillae):2:(L+L_ancillae))...)
        new(pbc, L, L_ancillae, oddfirst, sites, oddsites, evensites)
    end
end

# State container
mutable struct State{L}
    const L::Int64
    const indices::NTuple{L, Int64}
    const flatdim::Int64
    const tensordim::NTuple{L, Int64}
    state::Vector{ComplexF64}
    buffer::Vector{ComplexF64} # Reduces extra allocations in intermediate steps 
    perm::Vector{Int64} # Stores the permutation of sites during dynamics
    function State(L::Int64, state)
        perm = collect(1:L)
        new{L}(L, tuple(perm...), 2^L, tuple(fill(2, L)...), 
               state, copy(state), perm)
    end
end


# Swap array elements
function swap!(array, index1, index2)
    array[index1], array[index2] = array[index2], array[index1]
end

# Find right permutation in the permuted indices basis. 
function find_perm(sites, psi::State)
    Lgate = size(sites, 1)
    perm::Vector{Int64} = collect(psi.indices)
    state_perm = copy(psi.perm)
    for (i,site) in enumerate(sites)
        # Find the correct sites in the permuted indices basis
        perm_site = findfirst(isequal(site), state_perm)
        swap!(perm, i, perm_site)
        swap!(state_perm, i, perm_site)
    end
    return perm
end

function perm_sites!(psi::State, perm::Vector{Int64})
    permute!(psi.perm, perm)
    permutedims!(reshape(psi.buffer, psi.tensordim), 
                 reshape(psi.state, psi.tensordim), perm)
    psi.state, psi.buffer = psi.buffer, psi.state
end

function group_sites!(psi::State, sites::Vector{Int64})
    nsites = size(sites, 1)
    # Flatten in reverse order due to Julia is column majored arrays
    perm = find_perm(reverse(sites), psi) 
    permute!(psi.perm, perm)
    permutedims!(reshape(psi.buffer, psi.tensordim), 
                 reshape(psi.state, psi.tensordim), perm)
    return reshape(psi.buffer, (2^nsites, 2^(psi.L-nsites)))
end

# Applies any k-qubit gate to any k sites
function apply_gate!(psi::State, sites::Vector{Int64}, gate::Matrix{ComplexF64})
    permstate::Matrix{ComplexF64} = group_sites!(psi, sites)
    mul!(reshape(psi.state, size(permstate)), gate, permstate)
end

# Measurement operator container
struct MeasurementOperator
    Lop::Int64 
    operator::Hermitian{ComplexF64, Matrix{ComplexF64}}
    eigenvalues::Vector{Float64} 
    eigenbras::Vector{Matrix{ComplexF64}} 
    eigenkets::Vector{Matrix{ComplexF64}} 
    function MeasurementOperator(operator)
        op = Hermitian(ComplexF64.(operator))
        L = trailing_zeros(size(op, 1))

        decomp = eigen(op)
        vals = unique(x -> round(x, digits = 10), decomp.values)

        eigenbras = Vector{Matrix{ComplexF64}}()
        eigenkets = Vector{Matrix{ComplexF64}}()

        for (i,λ) in enumerate(vals)
            indices = findall(x -> isapprox(x, λ, atol=1E-9), decomp.values)
            push!(eigenkets, decomp.vectors[:,indices])
            push!(eigenbras, Matrix(eigenkets[end]'))
        end
        new(L, op, vals, eigenbras, eigenkets)
    end
end

## Generalized many-qubit projective measurement 
function measure!(psi::State, sites::Vector{Int64}, op::MeasurementOperator)
    permstate::Matrix{ComplexF64} = group_sites!(psi, sites)
    r::Float64, cumul::Float64 = rand(), 0
    for (i,bras) in enumerate(op.eigenbras)
        dim = (size(bras, 1), size(permstate, 2))
        mul!(reshape(view(psi.state, 1:prod(dim)), dim), bras, permstate)
        prob = norm(view(psi.state, 1:prod(dim)))^2
        cumul += prob

        if(cumul < r && i < length(op.eigenbras)) continue end
        mul!(reshape(psi.buffer, size(permstate)), op.eigenkets[i]/sqrt(prob), 
             reshape(view(psi.state, 1:prod(dim)), dim))
        psi.state, psi.buffer = psi.buffer, psi.state
        return (prob, op.eigenvalues[i])
    end

    #postmeas_states = [bras * permstate for bras in op.eigenbras]
    #probs = [norm(state)^2 for state in postmeas_states]
    #index = findfirst(>(rand() * sum(probs)), cumsum(probs))
    #psi.state = reshape((op.eigenkets[index] * postmeas_states[index]) / sqrt(probs[index]), 2^psi.L)
    #return (probs[index], op.eigenvalues[index])
end

# Localilly measures a qubit in the Z direction
function measure!(psi::State, site::Int64)
    perm_site = findfirst(isequal(site), psi.perm)
    state = reshape(psi.state, (2^(perm_site-1), 2, 2^(psi.L-perm_site)))
    prob = norm(@views state[:,1,:])^2
    if(rand() < prob)
        state[:,2,:] .= 0
        state ./= sqrt(prob)  
        res = (prob, 1)
    else
        state[:,1,:] .= 0
        state ./= sqrt(1-prob) 
        res = (1-prob, 0)
    end
    return res
end

# Applying gates might permute the site basis of the state
# This function sorts the sites to the original order
# Careful, not tested!
function sort_sites!(psi::State)
    perm_sites!(sortperm(psi.perm))
end

function overlap(bra::State, ket::State)
    bra_state = bra.state
    if(bra.perm != ket.perm)
        perm = find_perm(ket.perm, bra)
        permutedims!(reshape(bra.buffer, bra.tensordim), 
                     reshape(bra_state, bra.tensordim), perm) 
        bra_state = bra.buffer 
    end
    return dot(bra_state, ket.state)
end

# A usefull progress printer with timestamps
include("progress.jl")

# Predifined gates
include("gates.jl")

# Simple and useful many-body states
include("states.jl")

# Moments and frame potentials of the Haar ensemble
include("haar_moments.jl")

# Moments, frame potentials, and Haar distances of projected ensembles
include("projected_ensemble.jl")

# Vectors of sites for applications of gates in a Trotter decomposition
include("trotter_sites.jl")

# Observables
include("observables.jl")

# Apply gates and operations iteratively many times
include("iterative_applications.jl")

# module end
end
