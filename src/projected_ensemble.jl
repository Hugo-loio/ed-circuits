#
# Functions is this file compute ensembles of states projected in a complement of a subsystem A
# Further more they also compute moments of the ensembles, frame potentials and distances to the Haar ensemble
#

function exact_projected_ensemble(psi::State, sitesA::Vector{Int64})
    LA = size(sitesA, 1)
    LB::Int64 = psi.L - LA
    reshaped_state = reshape(psi.state, psi.tensordim)
    perm = find_perm(sitesA, psi)
    reshaped_state = reshape(permutedims(reshaped_state, perm), (2^LA, 2^LB))
    probs = Vector{Float64}([])
    states = Vector{Vector{ComplexF64}}([])
    for i in 1:2^LB
        prob = norm(reshaped_state[:,i])^2
        if(prob > 1E-13)
            push!(probs, prob)
            push!(states, reshaped_state[:, i] ./ sqrt(prob))
        end
    end
    return probs, states
end

# Group identical (up to a global phase) projected states 
function grouped_exact_projected_ensemble(psi::State, sitesA::Vector{Int64})
    LA = size(sitesA, 1)
    LB::Int64 = psi.L - LA
    reshaped_state = reshape(psi.state, psi.tensordim)
    perm = find_perm(sitesA, psi)
    reshaped_state = reshape(permutedims(reshaped_state, perm), (2^LA, 2^LB))
    probs = Vector{Float64}([])
    states = Vector{Vector{ComplexF64}}([])
    for i in 1:2^LB
        prob = norm(reshaped_state[:,i])^2
        if(prob < 1E-10)
            continue
        end
        state = reshaped_state[:,i] ./ sqrt(prob)
        save = true
        for (j,saved_state) in enumerate(states)
            overlap = abs2(dot(saved_state, state))
            if(abs(overlap - 1) < 1E-9)
                probs[j] += prob
                save = false
                break
            end
        end
        if(save)
            push!(probs, prob)
            push!(states, reshaped_state[:, i] ./ sqrt(prob))
        end 
    end
    return probs, states
end

# Same thing as the previous but also returns the overlaps
function grouped_exact_projected_ensemble2(psi::State, sitesA::Vector{Int64})
    LA = size(sitesA, 1)
    LB = psi.L - LA
    perm = find_perm(sitesA, psi)
    fullstate = reshape(permutedims(reshape(psi.state, psi.tensordim), perm), (2^LA, 2^LB))

    states = Vector{Vector{ComplexF64}}()
    overlaps = Vector{Vector{Float64}}() 
    probabilities = Vector{Float64}()

    for i in 1:2^LB
        save = true
        prob = norm(fullstate[:,i])^2
        if(prob < 1E-13)
            continue
        end
        state = fullstate[:,i] ./ sqrt(prob)
        push!(overlaps, Vector{Float64}())
        for (j,saved_state) in enumerate(states)
            overlap = abs2(dot(saved_state, state))
            if(abs(overlap - 1) < 1E-9)
                probabilities[j] += prob
                pop!(overlaps)
                save = false
                break
            end
            push!(overlaps[end], overlap)
        end
        if(save)
            push!(states, state)
            push!(probabilities, prob)
        end
    end

    return probabilities, states, overlaps
end

function projected_state(psi::State, sitesA::Vector{Int64}, bit_str::Int64)
    LA = size(sitesA, 1)
    LB = psi.L - LA
    reshaped_state = reshape(psi.state, psi.tensordim)
    perm = find_perm(sitesA, psi)
    reshaped_state = reshape(permutedims(reshaped_state, perm), (2^LA, 2^LB))
    projected_state = copy(reshaped_state[:,bit_str])
    prob = norm(projected_state)^2
    if(prob == 0)
        return prob, projected_state 
    end
    return prob, projected_state ./ sqrt(prob)
end

function random_projected_state(projected_state::Array{ComplexF64}, LA::Int64, LB::Int64)
    rands = rand(LB) 
    total_prob::Float64 = 1
    for i in 1:LB
        projected_state = reshape(projected_state, (2^(LA + LB - i), 2))
        prob = norm(projected_state[:,1])^2
        if rands[i] < prob
            projected_state = projected_state[:,1] ./ sqrt(prob)
            total_prob *= prob
        else
            projected_state = projected_state[:,2] ./ sqrt(1-prob)
            total_prob *= (1-prob)
        end
    end
    return total_prob, reshape(projected_state, 2^LA)
end

function random_projected_ensemble(psi::State, sitesA::Vector{Int64}, n::Int64)
    LA = size(sitesA, 1)
    LB = psi.L - LA
    perm = find_perm(sitesA, psi)
    permuted_state = permutedims(reshape(psi.state, psi.tensordim), perm)

    probs = zeros(Float64, n)
    states = Vector{Vector{ComplexF64}}()
    for i in 1:n
        probs[i], state = random_projected_state(permuted_state, LA, LB)
        push!(states, state)
    end
    return probs, states
end

function ensemble_moments(kmax::Int64, probs::Vector{Float64}, states::Vector{Vector{ComplexF64}}; moments::Vector{Array{ComplexF64,2}} = Vector{Array{ComplexF64,2}}([]))
    n = size(probs, 1)
    d = size(states[1], 1)
    if(length(moments) == 0)
        moments = [zeros(ComplexF64, (d^k, d^k)) for k in 1:kmax]
    end
    for i in 1:n
        rho = states[i] * adjoint(states[i])
        rhok = probs[i] * rho
        moments[1] .+= rhok
        for k in 2:kmax
            rhok = kron(rho, rhok)
            moments[k] .+= rhok
        end
    end
    return moments
end

function ensemble_moments(kmax::Int64, states::Vector{Vector{ComplexF64}}; moments::Vector{Array{ComplexF64,2}} = Vector{Array{ComplexF64,2}}([]))
    n = length(states)
    d = length(states[1])
    if(length(moments) == 0)
        moments = [zeros(ComplexF64, (d^k, d^k)) for k in 1:kmax]
    end
    for i in 1:n
        rho = states[i] * adjoint(states[i])
        rhok = rho / n
        moments[1] .+= rhok
        for k in 2:kmax
            rhok = kron(rho, rhok)
            moments[k] .+= rhok
        end
    end
    return moments
end

function mynorm(matrix::Matrix, normtype::String="frob")
    if (normtype == "frob")
        return norm(matrix)
    elseif (normtype == "trace")
        return sum(svdvals(matrix)) 
    else
        throw(ArgumentError(norm, "Norm not recognized"))
    end
end

function haar_dist(moments::Vector{Array{ComplexF64,2}}; normtype::String = "frob")
    ks = collect(1:size(moments,1))
    d = size(moments[1],1)
    return [mynorm(moments[k] - haar_moment(d,k), normtype) for k in ks]
end

function frame_potential(moments::Vector{Array{ComplexF64,2}})
    ks = collect(1:size(moments,1))
    return [norm(moments[k])^2 for k in ks]
end

# Haar distance from the moments of a random monte carlo projected ensemble
function rand_proj_moments_haar_dists(
        psi::State, sitesA::Vector{Int64}, kmax::Int64, nproj::Vector{Int}; 
        normtype::String = "frob",
        states::Vector{Vector{ComplexF64}} = Vector{Vector{ComplexF64}}([])
    )
    if(length(states) == 0)
        _, states = random_projected_ensemble(psi, sitesA, nproj[1])
    end
    moments = ensemble_moments(kmax, states)
    dists = [haar_dist(moments; normtype = normtype)]
    delta_n_old = nproj[1]
    for (i,n) in enumerate(nproj[2:end])
        delta_n = n - nproj[i]
        _, states = random_projected_ensemble(psi, sitesA, delta_n)
        moments = ensemble_moments(kmax, states; moments = moments*(delta_n_old/delta_n))
        push!(dists, haar_dist(moments*(delta_n/n); normtype = normtype))
        delta_n_old = delta_n
    end
    return dists
end

# Frame potential from the moments of a random Monte Carlo projected ensemble
function rand_proj_moments_frame_potential(
        psi::State, sitesA::Vector{Int64}, kmax::Int64, nproj::Vector{Int}; 
        normtype::String = "frob", 
        states::Vector{Vector{ComplexF64}} = Vector{Vector{ComplexF64}}([])
    )
    if(length(states) == 0)
        _, states = random_projected_ensemble(psi, sitesA, nproj[end])
    end
    moments = ensemble_moments(kmax, states[1:nproj[1]])
    frames = [[mynorm(moment, normtype)^2 for moment in moments]]
    delta_n_old = nproj[1]
    for (i,n) in enumerate(nproj[2:end])
        delta_n = n - nproj[i]
        moments = ensemble_moments(kmax, states[nproj[i]:n]; moments = moments*(delta_n_old/delta_n))
        push!(frames, [mynorm(moment*(delta_n/n), normtype)^2 for moment in moments])
        delta_n_old = delta_n
    end
    return frames
end

# Frame pontential from the overlaps of an exact projected ensemble
function exact_proj_overlaps_frame_potential(
        psi::State, sitesA::Vector{Int64}, kmax::Int64; 
        #normtype::String = "frob"
    )
    probs, states = exact_projected_ensemble(psi, sitesA)
    nstates = length(probs)
    frames = zeros(Float64, kmax)
    for i in 1:nstates
        frames .+= probs[i]*probs[i]
        # Off diagonal overlaps
        for j in i+1:nstates
            overlap = abs2(dot(states[i], states[j]))
            prob = probs[i]*probs[j]
            for k in 1:kmax
                frames[k] += 2*prob*overlap^k
            end
        end
    end
    return frames
end

# Frame potential of a state ensemble from the overlap formula
function overlaps_frame_potential(states::Vector{Vector{ComplexF64}}, kmax::Int64, ns::Vector{Int64})
    frames = Vector{Vector{Float64}}([])
    nstart = 1
    overlaps = zeros(Float64, kmax)
    ks = collect(1:kmax)
    for n in ns
        # Diagonal overlaps
        for i in nstart:n
            overlaps .+= 1
        end
        # Off-diagonal elements : top block
        for i in 1:nstart-1
            for j in nstart:n
                overlap = abs2(dot(states[i], states[j]))
                overlaps += 2*(overlap .^ ks)
            end
        end
        # Off-diagonal elements : diagonal block
        for i in nstart:n
            for j in i+1:n
                overlap = abs2(dot(states[i], states[j]))
                overlaps += 2*(overlap .^ ks)
            end
        end
        push!(frames, overlaps/n^2)
        nstart = n+1
    end
    return frames
end

# Frame potential from the overlaps of a random Monte Carlo projected ensemble
function rand_proj_overlaps_frame_potential(
        psi::State, sitesA::Vector{Int64}, kmax::Int64, nproj::Vector{Int64}; 
        states::Vector{Vector{ComplexF64}} = Vector{Vector{ComplexF64}}([])
    )
    if(length(states) == 0)
        _, states = random_projected_ensemble(psi, sitesA, nproj[end])
    end
    return overlaps_frame_potential(states, kmax, nproj)
end

# Only works if the projected ensemble measure is invariant under spin flips
function rand_proj_overlaps_frame_potential_v2(
        psi::State, sitesA::Vector{Int64}, kmax::Int64, nproj::Vector{Int}; 
        #normtype::String = "frob", 
        states::Vector{Vector{ComplexF64}} = Vector{Vector{ComplexF64}}([])
    )
    if(length(states) == 0)
        _, states = random_projected_ensemble(psi, sitesA, nproj[end])
    end

    LA = size(sitesA, 1)
    LB = psi.L - LA
    p0, psi0 = projected_state(psi, sitesA, 2^LB)
    prefactor = p0 * 2^LB

    frames = Vector{Vector{Float64}}([])
    nstart = 1
    overlaps = zeros(Float64, kmax)
    for n in nproj
        for i in nstart:n
            overlap = abs2(dot(psi0, states[i]))
            for k in 1:kmax
                overlaps[k] += overlap^k
            end
        end
        push!(frames, prefactor*overlaps/n)
        nstart = n
    end
    return frames
end

# Frame potential from overlaps and state multiplicity
# overlaps[i][j] gives overlap between state i and j, with i > j
function frame_potential_from_overlaps_multiplicity(
        overlaps::Vector{Vector{Float64}}, multiplicities::Vector{Int64}, 
        kmax::Int64, n::Int64
    )
    frames = zeros(Float64, kmax)
    ks = collect(1:kmax)
    for (i,m) in enumerate(multiplicities)
        # Identical state overlaps
        frames .+= m^2
        # Different state overlaps
        for (j,overlap) in enumerate(overlaps[i])
            frames += (2*m*multiplicities[j])*(overlap .^ ks)
        end
    end
    return frames/n^2
end

# Frame potential of a random projected ensemble, grouping identical states
function rand_proj_group_overlaps_frame_potential(
        psi::State, sitesA::Vector{Int64}, kmax::Int64, nproj::Vector{Int64}
    )
    LA = size(sitesA, 1)
    LB = psi.L - LA
    perm = find_perm(sitesA, psi)
    permuted_state = permutedims(reshape(psi.state, psi.tensordim), perm)

    states = Vector{Vector{ComplexF64}}()
    overlaps = Vector{Vector{Float64}}() 
    multiplicities = Vector{Int64}()
    frames = Vector{Vector{Float64}}()

    for n in 1:nproj[end]
        save = true
        _, state = random_projected_state(permuted_state, LA, LB)
        push!(overlaps, Vector{Float64}())
        for (j,saved_state) in enumerate(states)
            overlap = abs2(dot(saved_state, state))
            if(abs(overlap - 1) < 1E-9)
                multiplicities[j] += 1
                pop!(overlaps)
                save = false
                break
            end
            push!(overlaps[end], overlap)
        end
        if(save)
            push!(states, state)
            push!(multiplicities, 1)
        end

        if(n in nproj)
            push!(frames, frame_potential_from_overlaps_multiplicity(overlaps, multiplicities, kmax, n))
        end
    end

    return frames
end

# Frame potential from overlaps and probabilities 
# overlaps[i][j] gives overlap between state i and j, with i > j
function frame_potential_from_overlaps_probs(
        overlaps::Vector{Vector{Float64}}, probs::Vector{Float64}, kmax::Int64
    )
    frames = zeros(Float64, kmax)
    ks = collect(1:kmax)
    for (i,p) in enumerate(probs)
        # Identical state overlaps
        frames .+= p*p
        # Different state overlaps
        for (j,overlap) in enumerate(overlaps[i])
            frames += (2*p*probs[j])*(overlap .^ ks)
        end
    end
    return frames
end

# Frame potential of an exact projected ensemble, grouping identical states
function exact_proj_group_overlaps_frame_potential(
        psi::State, sitesA::Vector{Int64}, kmax::Int64
    )
    probabilities, _, overlaps = grouped_exact_projected_ensemble2(psi, sitesA)
    return frame_potential_from_overlaps_probs(overlaps, probabilities, kmax)
end

function safe_sqrt(x, verbose::Bool)
    if(x < 0 && verbose)
        @warn "Found a negative number in the sqrt, outputing sqrt(abs)... " x
    end
    return sqrt(abs(x))
end

function frame_potential_haar_dist(d::Int64, frames::Vector{Float64}; verbose::Bool=false)
    kmax = length(frames)
    haar_frames = [haar_frame_potential(d,k) for k in 1:kmax]
    [safe_sqrt(frames[k] - haar_frames[k], verbose) for k in 1:kmax] 
end

function frame_potential_haar_dist(d::Int64, frames::Vector{Vector{Float64}}; verbose::Bool=false)
    kmax = length(frames[1])
    haar_frames = [haar_frame_potential(d,k) for k in 1:kmax]
    [[safe_sqrt(frame[k] - haar_frames[k], verbose) for k in 1:kmax] for frame in frames]
end

