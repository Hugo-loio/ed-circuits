#
# This functions compute moments and frame potentials of an ensemble of random Haar states
#

function permutation_matrix(perm, d, k)
    dim = d^k
    mat = zeros(Int, dim, dim)
    for i in 0:(dim - 1)
        indices = reverse(digits(i, base=d, pad=k))
        permuted_indices = [indices[p] for p in perm]
        j = foldl((x, y) -> x * d + y, permuted_indices)
        mat[i+1, j+1] = 1
    end
    return mat
end

function haar_moment(d, k)
    if d == 1
        return 1
    end
    perms = collect(permutations(1:k))
    rho = zeros(ComplexF64, d^k, d^k)
    for perm in perms
        rho += permutation_matrix(perm, d, k)
    end
    return rho * prod(d:(k + d - 1))
end

function haar_frame_potential(d, k)
    return 1 / binomial(k + d - 1, k)
end
