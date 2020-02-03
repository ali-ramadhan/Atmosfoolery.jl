abstract type AbstractAcousticTimeStepper end

struct ExplicitAcousticTimeStepper <: AbstractAcousticTimeStepper end
struct VerticallyImplicitAcousticTimeStepper <: AbstractAcousticTimeStepper end

function acoustic_time_steps(rk3_iter, nₛ, Δt)
    rk3_iter == 1 && return 1,         Δt/3
    rk3_iter == 2 && return Int(nₛ/2), Δt/nₛ
    rk3_iter == 3 && return nₛ,        Δt/nₛ
end

AU() =

function acoustic_time_stepping!(grid)
    Nτ, Δτ = acoustic_time_steps(rk3_iter, nₛ, Δt)

    for n in 1:Nτ
        for k in 1:Nz, j in 1:Ny, i in 1:Nx
            @inbounds begin
                U″[i, j, k] = U″[i, j, k] + Δτ * AU()
                V″[i, j, k] = V″[i, j, k] + Δτ * AV()
                W″[i, j, k] = W″[i, j, k] + Δτ * AW()
                Θ″[i, j, k] = Θ″[i, j, k] + Δτ * AΘ()
                ρ″[i, j, k] = ρ″[i, j, k] + Δτ * Aρ()
            end
        end
    end
end
