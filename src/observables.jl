# 
# Some functions to compute observables
# 

function reduced_density_matrix(psi::State, sitesA)
    LA = size(sitesA, 1)
    reshaped_state = reshape(psi.state, psi.tensordim)
    perm = find_perm(sitesA, psi)
    reshaped_state = reshape(permutedims(reshaped_state, perm), (2^LA, 2^(psi.L-LA)))
    rhoA = reshaped_state * transpose(conj(reshaped_state))
    return rhoA
end

## This is not very efficient
function von_neumann_ent(rho)
    eigs = 1E-9 .+ real.(eigvals(rho))
    return -sum(eigs .* log.(eigs))
end

function bip_ent_entropy(psi::State, sitesA::Vector)
    LA = size(sitesA, 1)
    state = reshape(psi.state, psi.tensordim)
    perm = find_perm(sitesA, psi)
    ss = svdvals(reshape(permutedims(state, perm), (2^LA, 2^(psi.L-LA))))
    ps = 1E-9 .+ ss .* ss
    return -sum(ps .* log.(ps))
end

# Bipartition starts at site
function bip_ent_entropy(psi::State, site::Integer)
    return bip_ent_entropy(psi, collect(1:site-1))
end

function tripartite_mutual_information(psi::State, sitesA, sitesB, sitesC)
    SA = bip_ent_entropy(psi, sitesA)
    SB = bip_ent_entropy(psi, sitesB)
    SC = bip_ent_entropy(psi, sitesC)
    SAUB = bip_ent_entropy(psi, vcat(sitesA, sitesB))
    SAUC = bip_ent_entropy(psi, vcat(sitesA, sitesC))
    SBUC = bip_ent_entropy(psi, vcat(sitesB, sitesC))
    SAUBUC = bip_ent_entropy(psi, vcat(sitesA, sitesB, sitesC))
    SA + SB + SC - SAUB - SAUC - SBUC + SAUBUC
end

function local_magnetization(psi::State, site)
    perm_site = findfirst(isequal(site), psi.perm)
    state = reshape(psi.state, (2^(perm_site-1), 2, 2^(psi.L-perm_site)))
    return 2*norm(state[:,1,:])^2 - 1
end

function magnetization(psi::State, sites::Vector{Int} = Vector{Int}([]))
    if(length(sites) == 0)
        sites = collect(1:psi.L)
    end
    res = 0
    for site in sites
        res += local_magnetization(psi, site)
    end
    return res
end

# Rényi stabilizer entropy exact - O(4^L)
function stabilizer_entropy(psi::State, ns::Vector{Int})
    I = sparse([1 0; 0 1])
    X = sparse([0 1; 1 0])
    Y = sparse([0 -im; im 0])
    Z = sparse([1 0; 0 -1])
    paulis = [I,X,Y,Z]
    D = 2^psi.L

    Xis = zeros(Float64, 4^psi.L)
    for (i, string) in enumerate(Iterators.product(fill(1:4, psi.L)...))
        paulistring = reduce((a, b) -> kron(a, b), (paulis[j] for j in string)) 
        psi_p = paulistring * psi.state
        Xis[i] = real(dot(psi.state, psi_p))^2 / D
    end

    res = Vector{Float64}([])

    for n in ns
        if(n == 1)
            push!(res, -sum(Xis .* log.(Xis .+ 1e-10)) - psi.L * log(2))
        elseif (n > 1)
            push!(res, log(sum(Xis.^n)) / (1 - n) - psi.L * log(2))
        end
    end

    return res
end

# Expectation value of a Pauli string over sites
function pauli_string_ev(psi::State, string::String, sites::Vector{Int64}) 
    flip_mask = phase_mask = num_y = 0
    
    perm_sites = [findfirst(isequal(site), psi.perm) for site in sites]
    for (i, s) in enumerate(perm_sites)
        pauli = lowercase(string[i])
        mask = 1 << (s - 1)
        if pauli == 'x' 
            flip_mask |= mask
        elseif pauli == 'y' # Y = iXZ
            flip_mask |= mask
            phase_mask |= mask
            num_y += 1
        elseif pauli == 'z'
            phase_mask |= mask
        end
    end
    
    total_ev = 0.0 + 0.0im

    @inbounds for i in 0:psi.flatdim-1
        target = i ⊻ flip_mask # XOR is a controled not (flips the right indices) 
        phase = iseven(count_ones(target & phase_mask)) ? 1.0 : -1.0
        total_ev += conj(psi.state[i+1]) * phase * psi.state[target+1]
    end
    return real(total_ev * (1im)^num_y)
end
