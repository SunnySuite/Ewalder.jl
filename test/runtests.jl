using TestItemRunner

@run_package_tests


@testitem "Madelung Constants" begin
    using LinearAlgebra

    # Madelung constants of ionic crystals:
    # - Y. Sakamoto, The J. of Chem. Phys., vol. 28, (1958), p. 164
    E_NaCl = -1.7475645946331822
    E_CsCl = -1.76267477307099

    # CsCl
    latvecs = [[1,0,0], [0,1,0], [0,0,1]]
    pos = [[0,0,0], [0.5,0.5,0.5]]
    sys = Ewalder.System(; latvecs, pos)
    E = Ewalder.energy(sys; charges=[1., -1.])
    r = norm(pos[2]-pos[1])
    @test isapprox(E*r, E_CsCl; atol=1e-13)

    # CsCl, upscaled by 2
    sys = Ewalder.System(; latvecs=2latvecs, pos=2pos)
    E = Ewalder.energy(sys; charges=[1., -1.])
    @test isapprox(E*2r, E_CsCl; atol=1e-13)

    # CsCl, tilted cell
    latvecs = [[1,0,0], [1,1,0], [0,0,1]]
    pos = [[0,0,0], [0.5,0.5,0.5]]
    sys = Ewalder.System(; latvecs, pos)
    E = Ewalder.energy(sys; charges=[1., -1.])
    r = norm(pos[2]-pos[1])
    @test isapprox(E*r, E_CsCl; atol=1e-13)

    # NaCl, cubic unit cell
    latvecs = [[2,0,0], [0,2,0], [0,0,2]]
    pos = Tuple.(CartesianIndices((0:1,0:1,0:1)))[:]
    parity(x) = 2mod(x, 2) - 1
    charges = [parity(p[1]+p[2]+p[3]) for p=pos]
    sys = Ewalder.System(; latvecs, pos)
    E = Ewalder.energy(sys; charges)
    @test isapprox(E/4, E_NaCl; atol=1e-13)
    
    # NaCl, primitive cell
    latvecs = [[1,1,0], [1,0,1], [0,1,1]]
    pos = [[0,0,0], [1,1,1]]
    charges = [1, -1]
    sys = Ewalder.System(; latvecs, pos)
    E = Ewalder.energy(sys; charges)
    @test isapprox(E, E_NaCl; atol=1e-13)

end


# Accuracy should be fairly consistent, given a fixed c0 parameter.
@testitem "Accuracy consistency" begin
    using LinearAlgebra

    E_CsCl = -1.76267477307099
    c0 = 3

    pos = [[0,0,0], [0.5,0.5,0.5]]
    r = norm(pos[2]-pos[1])
    
    # Reference calculation
    latvecs = [[1,0,0], [0,1,0], [0,0,1]]
    sys = Ewalder.System(; latvecs, pos, c0, c1=2)
    ΔE1 = E_CsCl - Ewalder.energy(sys; charges=[1., -1.])*r
    
    # Sheared box
    latvecs = [[1,0,0], [10,1,0], [0,0,1]]
    sys = Ewalder.System(; latvecs, pos, c0, c1=2)
    ΔE2 = E_CsCl - Ewalder.energy(sys; charges=[1., -1.])*r
    
    # Large Fourier space weight
    latvecs = [[1,0,0], [0,1,0], [0,0,1]]
    sys = Ewalder.System(; latvecs, pos, c0, c1=0.4)
    ΔE3 = E_CsCl - Ewalder.energy(sys; charges=[1., -1.])*r
    
    # Large real space weight
    latvecs = [[1,0,0], [0,1,0], [0,0,1]]
    sys = Ewalder.System(; latvecs, pos, c0, c1=10)
    ΔE4 = E_CsCl - Ewalder.energy(sys; charges=[1., -1.])*r
    
    # For the parameters chosen above, all errors are on the same scale. This test
    # is not very robust, and might fail after internal changes to the Ewald sum. In
    # this case, some hand-testing could be required to pick new parameters.
    @test all([ΔE1, ΔE2, ΔE3, ΔE4]) do x
        1e-5 < abs(x) < 1e-3
    end
end


# For a system with both charges and dipoles, the Ewald result should be largely
# independent of the parameter c1 which balances between real and Fourier space
# summation.
@testitem "Real/Fourier space consistency" begin
    using LinearAlgebra

    latvecs = [[1,0,0], [0,1,0], [0,0,1]]
    pos = [[0,0,0], [0.6,0.4,0.3]] # broken symmetry
    charges = [1., -1.]
    dipoles = [[0.35, -0.27, 0.8], [-0.1, 0.5, 0.32]] # arbitrary choice

    sys = Ewalder.System(; latvecs, pos, c1=2)
    E_ref = Ewalder.energy(sys; charges, dipoles)

    sys = Ewalder.System(; latvecs, pos, c1=0.4)
    E = Ewalder.energy(sys; charges, dipoles)
    @test isapprox(E, E_ref; atol=1e-11)

    sys = Ewalder.System(; latvecs, pos, c1=10)
    E = Ewalder.energy(sys; charges, dipoles)
    @test isapprox(E, E_ref; atol=1e-11)
end


# A dipole should be effectively the same as two slightly displaced,
# opposite-signed charges.
@testitem "Dipoles as finite differences" begin
    using LinearAlgebra

    ### Reference energy of a single periodic dipole
    
    latvecs = [[1,0,0], [0,1,0], [0,0,1]]
    p = [0.2, 0.5, 0.7]
    r = [0.5, 0.5, 0.5]
    sys = Ewalder.System(; latvecs, pos=[r])
    E_ref = Ewalder.energy(sys; dipoles=[p])
    
    ### Perform same calculation using finite differences
    
    ϵ = 0.005 # carefully tuned to minimize FP roundoff
    r1 = r - ϵ*p/2
    r2 = r + ϵ*p/2
    charges = [1, -1] * (1/ϵ)
    sys = Ewalder.System(; latvecs, pos=[r1, r2])
    E_ewald = Ewalder.energy(sys; charges)
    
    Δr = ϵ*norm(p)
    E_internal = charges[1]*charges[2]/Δr
    
    E_fd = E_ewald - E_internal
    
    ### Check answer to 4 digits. Cannot expect too much accuracy here.

    @test isapprox(E_ref, E_fd; atol=1e-4)    
end
