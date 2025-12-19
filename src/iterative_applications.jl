# 
# Functions to apply gates many times or in layers
# 

function measure_layer!(circuit::Circuit, psi::State, prob::Float64)
    rands = rand(Float64, circuit.L)
    sites = circuit.sites[findall(x -> x < prob, rands)]
    record = []
    for site in sites
        push!(record, measure!(psi, site))
    end
    return record
end

function apply_gate_manytimes!(psi::State, sites_vector::Vector{Vector{Int64}}, gate::Matrix{ComplexF64})
    for sites in sites_vector
        apply_gate!(psi, sites, gate)
    end
end

function apply_gates!(psi::State, sites_vector::Vector{Vector{Int64}}, gates::Vector{Matrix{ComplexF64}})
    for (i,sites) in enumerate(sites_vector)
        apply_gate!(psi, sites, gates[i])
    end
end

function apply_random_gates!(psi::State, sites_vector::Vector{Vector{Int64}}, gates::Vector{Matrix{ComplexF64}})
    n = size(gates,1)
    for sites in sites_vector
        apply_gate!(psi, sites, gates[rand(1:n)])
    end
end

function apply_random_haar_gates!(psi::State, sites_vector::Vector{Vector{Int64}})
    for sites in sites_vector
        nsites = size(sites,1)
        apply_gate!(psi, sites, haar_gate(nsites))
    end
end
