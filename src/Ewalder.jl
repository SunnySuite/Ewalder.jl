module Ewalder

import StaticArrays: SVector, SMatrix, SA
import LinearAlgebra: normalize, ⋅, ×, inv
import SpecialFunctions: erfc

const Vec3 = SVector{3,Float64}
const Mat3 = SMatrix{3,3,Float64,9}

# Two ion indices i and j, and a cell offset n
Base.@kwdef struct Neighbor
    i::Int32
    j::Int32
    n::SVector{3,Int32}
end

"""
    System(; latvecs::Vector{NTuple{3, Float}},
             pos::Vector{Vec3}
             c0::Float64 = 6.0
             c1::Float64 = 2.0)

Create a system that is periodic in each of the three supercell lattice vectors
`latvecs` and with atom positions `pos`. The optional parameter `c0` controls
accuracy of Ewald summation. The default is `c0 = 6`, which yields errors of
order `1e-12` (smaller `c0` will run faster). The optional parameter `c1`
controls the balance between real and Fourier space summation (larger `c1`
implies more real-space summation).
"""
Base.@kwdef struct System
    latvecs::SVector{3, Vec3}
    pos::Vector{Vec3}
    c0::Float64 = 6.0 # Controls accuracy
    c1::Float64 = 2.0 # Controls balance of real vs Fourier space
end

function recipvecs(sys::System)
    b = 2π*inv(reduce(hcat, sys.latvecs))
    SA[b[1,:], b[2,:], b[3,:]]
end

function volume(sys::System)
    a = sys.latvecs
    abs((a[1]×a[2])⋅a[3])
end

function displacement(sys::System, neigh::Neighbor)
    (; i, j, n) = neigh
    ri = sys.pos[i]
    rj = sys.pos[j]
    return rj - ri + n'*sys.latvecs
end

function distance2(sys::System, neigh::Neighbor)
    r = displacement(sys, neigh)
    return r⋅r
end

function sigma(sys::System)
    N = length(sys.pos)
    N == 0 && error("Empty system")
    L = cbrt(volume(sys))
    return L / (sys.c1 * N^(1/6))
end

real_space_cutoff(sys::System)    = √2 * sys.c0 * sigma(sys)
fourier_space_cutoff(sys::System) = √2 * sys.c0 / sigma(sys)


function wrap_positions!(sys::System; warn=false)
    # Lattice vectors as a 3x3 matrix
    A = reduce(hcat, sys.latvecs)
    invA = inv(A)

    for (i,r) = enumerate(sys.pos)
        # Convert r to basis of lattice vectors
        v = invA * r

        # Optionally warn user
        inrange(x) = 0 ≤ x ≤ 1
        if warn && !all(inrange, v)
            println("Warning: Ion at $r is outside unit cell. Its fractional coordinates are $v.")
        end

        # Wrap coordinates
        sys.pos[i] = A * Vec3(mod(x, 1) for x = v)
    end
end

# Perform O(N^2) calculation of neighbor list.
function get_neighbors(sys::System; rmax=0.0)
    # Make sure positions are valid
    wrap_positions!(sys; warn=true)

    # Real space cutoff
    if iszero(rmax)
        rmax = real_space_cutoff(sys)
    end

    # Maximum latvec displacements, in each direction
    as = sys.latvecs
    bs = recipvecs(sys)
    nmax = map(as, bs) do a, b
        round(Int, rmax / (a⋅normalize(b)) + 1e-6) + 1
    end

    N = length(sys.pos)

    ret = Neighbor[]
    for i = 1:N, j = 1:N
        for n1 = -nmax[1]:nmax[1], n2 = -nmax[2]:nmax[2], n3 = -nmax[3]:nmax[3]
            n = (n1, n2, n3)
            neigh = Neighbor(; i, j, n)
            if distance2(sys, neigh) <= rmax*rmax
                if i != j || !iszero(Vec3(n))
                    push!(ret, neigh)
                end
            end
        end
    end
    return ret
end


@doc raw"""
    energy(sys::System; charges=Float64[], dipoles=Vec3[])

Calculate the Ewald energy in units of $1/4\pi\epsilon_0$ for a periodic system
of `charges` and `dipoles`. If either charge or dipole parameters are omitted,
they will be assumed zero.
"""
function energy(sys::System; neighbors=Neighbor[], charges=Float64[], dipoles=Vec3[])
    N = length(sys.pos)

    # Find neighbors if not provided
    if isempty(neighbors)
        neighbors = get_neighbors(sys)
    end

    # Set charges and dipoles to zero if not provided
    charges = isempty(charges) ? zeros(Float64, N) : convert(Vector{Float64}, charges)
    dipoles = isempty(dipoles) ? zeros(Vec3, N)    : convert(Vector{Vec3}, dipoles)

    # Check charge neutrality
    Q = sum(charges)
    @assert abs(Q) < 1e-12
    charges .-= Q / N

    # Energy accumulator
    E = 0.0

    #####################################################
    ## Real space sum
    σ = sigma(sys)
    σ² = σ*σ

    for neigh = neighbors
        (; i, j) = neigh
        
        qᵢ, pᵢ = (charges[i], dipoles[i])
        qⱼ, pⱼ = (charges[j], dipoles[j])
        rvec = displacement(sys, neigh)
        r² = rvec⋅rvec
        r = √r²
        r³ = r²*r
        rhat = rvec/r

        @assert abs(r) > 1e-12 "Detected zero-distance neighbor, $neigh."

        erfc0 = erfc(r/(√2*σ))
        gauss0 = √(2/π) * (r/σ) * exp(-r²/2σ²)

        # Charge-charge
        E += (1/2) * qᵢ * qⱼ * erfc0 / r

        # Charge-dipole
        E += (1/2) * (qᵢ*pⱼ - pᵢ*qⱼ) ⋅ (rhat/r²) * (erfc0 + gauss0)

        # Dipole-dipole
        E += (1/2) * ((pᵢ⋅pⱼ/r³) * (erfc0 + gauss0) - (3(pᵢ⋅rhat)*(pⱼ⋅rhat)/r³) * (erfc0 + (1+r²/3σ²) * gauss0))
    end

    #####################################################
    ## Fourier space sum
    V = volume(sys)
    as = sys.latvecs
    bs = recipvecs(sys)

    # Maximum recipvec displacements, in each direction
    kmax = fourier_space_cutoff(sys)
    mmax = map(as, bs) do a, b
        round(Int, kmax / (normalize(a)⋅b) + 1e-6) + 1
    end
    
    # Loop over grid of k points
    for m1 = -mmax[1]:mmax[1], m2 = -mmax[2]:mmax[2], m3 = -mmax[3]:mmax[3]
        k = Vec3(m1, m2, m3)' * bs
        k² = k⋅k
        if 0 < k² <= kmax*kmax
            ρk = 0.0im
            for (r, q, p) = zip(sys.pos, charges, dipoles)
                ρk += (q + im*(p⋅k)) * cis(-k⋅r)
            end
            E += 4π * (1/2V) * (exp(-σ²*k²/2) / k²) * abs2(ρk)
        end
    end


    #####################################################
    ## Remove self energies
    σ³ = σ^3
    for (q, p) = zip(charges, dipoles)
        E += - q*q/(√(2π)*σ) - (p⋅p)/(√(2π)*(3σ³))
    end

    return E
end


export System, energy


end
