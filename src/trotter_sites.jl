#
# Functions in this file produce a vector of sites where gates should act in a time step of a Trotter decomposition
# 


# First order Trotter decomposition
# For local any qubit gates
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
