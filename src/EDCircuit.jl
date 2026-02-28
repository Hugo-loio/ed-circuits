#
# This module implements exact diagonalization simulations for 1D spin-1/2 chains
#

module EDCircuit

using Dates
using LinearAlgebra
using SparseArrays
using SpecialFunctions: factorial
using Combinatorics: permutations

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

mutable struct State{L}
    const L::Int64
    const indices::NTuple{L, Int64}
    const flatdim::Int64
    const tensordim::NTuple{L, Int64}
    state::Vector{ComplexF64}
    perm::Vector{Int64} # Stores the permutation of sites during dynamics
    function State(L::Int64, state)
        perm = collect(1:L)
        new{L}(L, tuple(perm...), 2^L, tuple([2 for _ in 1:L]...), state, perm)
    end
end

function swap!(array, index1, index2)
    temp = array[index1]
    array[index1] = array[index2]
    array[index2] = temp
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

function apply_gate!(psi::State, sites::Vector{Int64}, gate::Matrix{ComplexF64})
    Lgate = size(sites, 1)
    perm = find_perm(reverse(sites), psi) # Since Julia is column majored the sites should be flattened in reverse order 
    # Permute sites and reshape state for multiplication
    psi.perm = psi.perm[perm]
    permstate::Matrix{ComplexF64} = reshape(permutedims(reshape(psi.state, psi.tensordim), perm), (2^Lgate, 2^(psi.L-Lgate)))
    #display(permstate)
    psi.state = reshape(gate * permstate, 2^psi.L)
end

function measure!(psi::State, site::Int64)
    perm_site = findfirst(isequal(site), psi.perm)
    state = reshape(psi.state, (2^(perm_site-1), 2, 2^(psi.L-perm_site)))
    prob = norm(state[:,1,:])^2
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
function sort_sites!(psi::State)
    tensor_state = reshape(psi.state, psi.tensordim)
    perm_state = permutedims(tensor_state, sortperm(psi.perm))
    psi.state = reshape(perm_state, 2^psi.L)
    psi.perm = collect(psi.indices)
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
