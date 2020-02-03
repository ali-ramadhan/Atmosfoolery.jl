using JULES.Operators

using Oceananigans: datatuple
using Oceananigans.BoundaryConditions: NoPenetrationBC, fill_halo_regions!
import Oceananigans: time_step!

####
#### Dirty hacks!
####

const hpbcs = HorizontallyPeriodicBCs()
const hpbcs_np = HorizontallyPeriodicBCs(top=NoPenetrationBC(), bottom=NoPenetrationBC())

####
#### Utilities for time stepping
####

function rk3_time_step(rk3_iter, Δt)
    rk3_iter == 1 && return Δt/3
    rk3_iter == 2 && return Δt/2
    rk3_iter == 3 && return Δt
end

####
#### Time-stepping algorithm
####

function time_step!(model::CompressibleModel; Δt, Nt=1)
    arch = model.architecture
    grid = model.grid
    time = model.clock.time
    coriolis = model.coriolis
    buoyancy = model.buoyancy
    closure = model.closure
    prog_temp = model.prognostic_temperature
    microphysics = model.microphysics
    forcing = model.forcing
    params = model.parameters

    Ũ  = model.momenta
    ρᵈ = model.density
    C̃  = model.tracers
    K̃ = model.diffusivities
    F  = model.slow_forcings
    R  = model.right_hand_sides
    IV = model.intermediate_vars

    pₛ = model.reference_pressure
    g  = model.gravity

    # On third RK3 step, we update Φ⁺ instead of model.intermediate_vars
    Φ⁺ = merge(Ũ, C̃, (ρ=ρᵈ,))

    # On the first and second RK3 steps we want to update intermediate Ũ and C̃.
    Ũ_names = propertynames(Ũ)
    IV_Ũ_vals = [getproperty(IV, U) for U in Ũ_names]
    IV_Ũ = NamedTuple{Ũ_names}(IV_Ũ_vals)

    C̃_names = propertynames(C̃)
    IV_C̃_vals = [getproperty(IV, C) for C in C̃_names]
    IV_C̃ = NamedTuple{C̃_names}(IV_C̃_vals)

    for _ in 1:Nt
        @debug "Computing slow forcings..."
        fill_halo_regions!(ρᵈ.data, hpbcs, arch, grid)
        fill_halo_regions!(datatuple(merge(Ũ, C̃)), hpbcs, arch, grid)
        fill_halo_regions!(Ũ.ρw.data, hpbcs_np, arch, grid)
        compute_slow_forcings!(F, grid, coriolis, closure, Ũ, ρᵈ, C̃, K̃, forcing, time, params)
        fill_halo_regions!(F.ρw.data, hpbcs_np, arch, grid)

        # RK3 time-stepping
        for rk3_iter in 1:3
            @debug "RK3 step #$rk3_iter..."

            @debug "  Computing right hand sides..."
            if rk3_iter == 1
                compute_rhs_args = (R, grid, prog_temp, buoyancy, microphysics, pₛ, g, ρᵈ, Ũ, C̃, F)
                fill_halo_regions!(ρᵈ.data, hpbcs, arch, grid)
                fill_halo_regions!(datatuple(merge(Ũ, C̃)), hpbcs, arch, grid)
                fill_halo_regions!(Ũ.ρw.data, hpbcs_np, arch, grid)
            else
                compute_rhs_args = (R, grid, prog_temp, buoyancy, microphysics, pₛ, g, IV.ρ, IV_Ũ, IV_C̃, F)
                fill_halo_regions!(IV.ρ.data, hpbcs, arch, grid)
                fill_halo_regions!(datatuple(merge(IV_Ũ, IV_C̃)), hpbcs, arch, grid)
                fill_halo_regions!(IV_Ũ.ρw.data, hpbcs_np, arch, grid)
            end

            compute_right_hand_sides!(compute_rhs_args...)

            # n, Δτ = acoustic_time_steps(rk3_iter)
            # acoustic_time_stepping!(Ũ, ρ, C, F, R; n=n, Δτ=Δτ)

            @debug "  Advancing variables..."
            LHS = rk3_iter == 3 ? Φ⁺ : IV
            advance_variables!(LHS, grid, Ũ, C̃, ρᵈ, R; Δt=rk3_time_step(rk3_iter, Δt))
        end

        model.clock.iteration += 1
        model.clock.time += Δt
    end

    return nothing
end
