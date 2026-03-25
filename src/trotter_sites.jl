#
# Functions in this file produce a vector of sites where gates should act in a time step of a Trotter decomposition
# 


# First order Trotter decomposition
# For local any qubit gates
# Deprecated, TODO: remove
function trotter_sites_vector(circuit::Circuit, sites_per_gate::Int64)
    sites_vector = Vector{Int64}[]
    index = 1
    layer = 1
    shift = sites_per_gate - 1
    while (layer <= sites_per_gate)
        if(circuit.sites[index] + shift > circuit.sites[end])
            if(circuit.pbc)
                sitesend = collect(circuit.sites[index:end])
                sitesbegin = collect(circuit.sites[1:index+shift-circuit.L])
                push!(sites_vector, vcat(sitesend, sitesbegin))
            end
            layer += 1
            index = layer
            continue
        end
        push!(sites_vector, collect(circuit.sites[index:index+shift]))
        if(circuit.sites[index] + shift == circuit.sites[end])
            layer += 1
            index = layer
        else
            index += sites_per_gate
        end
    end
    return sites_vector
end

# First order Trotter decomposition
# For long range 2-qubit gates (i, i + range)
# Deprecated, TODO: remove
function trotter_long_range_sites_vector(circuit::Circuit, range::Int64)
    sites_vector = Vector{Int64}[]
    index = 1
    layer = 1
    while (layer <= 2)
        if(index > circuit.L)
            index = 1 + layer*range
            layer += 1 
            continue
        end
        for _ in 1:range
            if(index > circuit.L)
                break
            end
            if(circuit.sites[index] + range > circuit.sites[end])
                if(circuit.pbc)
                    push!(sites_vector, [circuit.sites[index], circuit.sites[index+range-circuit.L]])
                end
            else
                push!(sites_vector, [circuit.sites[index], circuit.sites[index + range]])
            end
            index += 1 
        end
        if(index <= circuit.L)
            index += range
        end
    end
    return sites_vector
end

# Optimal gate site order for a first-order 1D Trotter decomposition
# for TI k-site interaction model ∑_i gate_(i,i + range, ..., i + (k-1)*range)
# Due to PBC, some edge gates in the same layer might not commute
function trotter_sites(L::Int64; 
        range::Int64 = 1, k::Int64 = 2, pbc::Bool = true)
    sites_vector = Vector{Vector{Int64}}()
    layers = Vector{UnitRange{Int64}}()
    layer_start = 1
    for layer in 1:k
        site = 1 + (layer-1)*range
        step = 1
        while(site <= L)
            sites = [mod1(site + x*range, L) for x in 0:k-1]
            if((site + (k-1)*range > L) && !pbc) break end
            push!(sites_vector, sites)
            site += (step % range == 0) ? (k - 1) * range + 1 : 1
            step += 1
        end
        push!(layers, layer_start:length(sites_vector))
        layer_start = length(sites_vector) + 1
    end
    return sites_vector, layers
end

# Generalization of the previous function for multiple dimensions
# assuming isotropic interactions in range and number of qubit
# This function applies a row wise flattening of the site indices
function trotter_sites(geometry::Vector{Int64}; 
        range::Int64 = 1, k::Int64 = 2, pbc::Bool = true)
    # Collection of groups of sites in lenght(geometry) dimensions
    sites_vector = Vector{Vector{Vector{Int64}}}() 
    layers = Vector{UnitRange{Int64}}()
    # Loop over geometry directions (gates in that direction)
    for dir in 1:length(geometry)
        partial_geometry = (geometry[1:dir-1]..., geometry[dir+1:end]...)
        sites_1d, layers_1d = trotter_sites(geometry[dir], range = range, 
                                            k = k, pbc = pbc)
        for layer_1d in layers_1d
            start_index = length(sites_vector) + 1
            for partial_site in CartesianIndices(partial_geometry)
                partial_site = collect(Tuple(partial_site))
                for sites in sites_1d[layer_1d]
                    push!(sites_vector, 
                          [vcat(partial_site[1:dir-1], site, 
                                partial_site[dir:end]) for site in sites])
                end
            end
            push!(layers, start_index:length(sites_vector))
        end
    end

    flat_rescale = ones(Int64, length(geometry))
    for i in length(geometry)-1:-1:1
        flat_rescale[i] = flat_rescale[i+1] * geometry[i+1]
    end
    flat_sites_vector = Vector{Vector{Int64}}() 
    for sites in sites_vector
        push!(flat_sites_vector, [sum((site .- 1) .* flat_rescale) + 1
                                  for site in sites])
    end

    return flat_sites_vector, layers
end
