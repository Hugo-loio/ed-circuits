#
# Optimized gate applications for pauli operators
#

# Local pauli gate
function apply_pauli!(pauli::Int64, site::Int64, psi::State)
    if pauli == 0 return end # Identity
    perm_site = findfirst(isequal(site), psi.perm)
    buffer = reshape(psi.buffer, (2^(perm_site-1), 2, 2^(psi.L - perm_site)))
    state = reshape(psi.state, (2^(perm_site-1), 2, 2^(psi.L - perm_site)))
    if pauli ==  1 # σ_x
        @views buffer[:, 1, :] .= state[:, 2, :]
        @views buffer[:, 2, :] .= state[:, 1, :]
        psi.state, psi.buffer = psi.buffer, psi.state
    elseif pauli == 2 # σ_y
        @views buffer[:, 1, :] .= -im .* state[:, 2, :]
        @views buffer[:, 2, :] .=  im .* state[:, 1, :]
        psi.state, psi.buffer = psi.buffer, psi.state
    elseif pauli == 3 # σ_z
        @views state[:,2,:] .*= -1
    end
end

# Pauli string gate
function apply_pauli_string!(string::Int64, sites::Vector{Int64}, psi::State)
    if(length(sites) > 31)
        error("Too many sites: string label can overflow")
    end
    for (i,site) in enumerate(sites)
        apply_pauli!(string % 4, site, psi)
        string >>= 2 # Bitwise shift by 2, identical to integer division by 4
    end
end

# Random Pauli string gate
function apply_random_pauli_string!(sites::Vector{Int64}, psi::State; 
        include_identity::Bool = false)
    min_str::Int64 = 1
    if(include_identity) min_str = 0 end
    max_str::Int64 = (1 << (2 * length(sites))) - 1
    apply_pauli_string!(rand(min_str:max_str), sites, psi)
end
