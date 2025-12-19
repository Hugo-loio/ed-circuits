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

# Bipartition starts at site
function bip_ent_entropy(psi::State, site::Integer)
    ss = svdvals(reshape(psi.state, (2^(site-1), 2^(psi.L-site+1))))
    ps = 1E-9 .+ ss .* ss
    return -sum(ps .* log.(ps))
end

function bip_ent_entropy(psi::State, sitesA::Vector)
    LA = size(sitesA, 1)
    state = reshape(psi.state, psi.tensordim)
    perm = find_perm(sitesA, psi)
    ss = svdvals(reshape(permutedims(state, perm), (2^LA, 2^(psi.L-LA))))
    ps = 1E-9 .+ ss .* ss
    return -sum(ps .* log.(ps))
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
