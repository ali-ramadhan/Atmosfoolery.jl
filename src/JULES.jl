module JULES

using OffsetArrays

T = Float64

#
# Grid definitions
#

struct Grid
    Nx
    Ny
    Nz
    Hx
    Hy
    Hz
    Δx
    Δy
    Δz
end

struct Field{X,Y,Z}
    data
    grid
end

#
# Field definitions
#

struct Cell end
struct Face end

function CellField(grid)
    data = OffsetArray{T}(
                          undef, 
                          -grid.Hx+1:grid.Nx+grid.Hx, 
                          -grid.Hy+1:grid.Ny+grid.Hy,
                          -grid.Hz+1:grid.Nz+grid.Hz
                         )
    Field{Cell,Cell,Cell}(data, grid)
end

function FaceFieldX(grid)
    data = OffsetArray{T}(
                          undef,
                          -grid.Hx+1:grid.Nx+grid.Hx,
                          -grid.Hy+1:grid.Ny+grid.Hy,
                          -grid.Hz+1:grid.Nz+grid.Hz
                         )
    Field{Face,Cell,Cell}(data, grid)
end

function FaceFieldY(grid)
    data = OffsetArray{T}(
                          undef,
                          -grid.Hx+1:grid.Nx+grid.Hx,
                          -grid.Hy+1:grid.Ny+grid.Hy,
                          -grid.Hz+1:grid.Nz+grid.Hz
                         )
    Field{Cell,Face,Cell}(data, grid)
end

function FaceFieldZ(grid)
    data = OffsetArray{T}(
                          undef,
                          -grid.Hx+1:grid.Nx+grid.Hx,
                          -grid.Hy+1:grid.Ny+grid.Hy,
                          -grid.Hz+1:grid.Nz+grid.Hz
                         )
    Field{Cell,Cell,Face}(data, grid)
end
  
@propagate_inbounds getindex(f::Field, inds...) = @inbounds getindex(f.data, inds...)
@propagate_inbounds setindex!(f::Field, v, inds...) = @inbounds setindex!(f.data, v, inds...)

@inline ∂xᶜᵃᵃ(i, j, k, grid, u) = @inbounds (u[i+1, j, k] - u[i,   j, k]) / grid.Δx
@inline ∂xᶠᵃᵃ(i, j, k, grid, c) = @inbounds (c[i,   j, k] - c[i-1, j, k]) / grid.Δx
@inline ∂yᵃᶜᵃ(i, j, k, grid, v) = @inbounds (v[i, j+1, k] - v[i, j,   k]) / grid.Δy
@inline ∂xᵃᶠᵃ(i, j, k, grid, c) = @inbounds (c[i, j,   k] - c[i, j-1, k]) / grid.Δy
@inline ∂xᵃᵃᶜ(i, j, k, grid, w) = @inbounds (w[i, j, k+1] - w[i, j, k  ]) / grid.Δz
@inline ∂xᵃᵃᶠ(i, j, k, grid, c) = @inbounds (c[i, j,   k] - c[i, j, k-1]) / grid.Δz

@inline Ixᶜᵃᵃ(i, j, k, grid, u) = @inbounds T(0.5) * (u[i+1, j, k] - u[i,   j, k]) 
@inline Ixᶠᵃᵃ(i, j, k, grid, c) = @inbounds T(0.5) * (c[i,   j, k] - c[i-1, j, k]) 
@inline Iyᵃᶜᵃ(i, j, k, grid, v) = @inbounds T(0.5) * (v[i, j+1, k] - v[i, j,   k]) 
@inline Ixᵃᶠᵃ(i, j, k, grid, c) = @inbounds T(0.5) * (c[i, j,   k] - c[i, j-1, k]) 
@inline Ixᵃᵃᶜ(i, j, k, grid, w) = @inbounds T(0.5) * (w[i, j, k+1] - w[i, j, k  ]) 
@inline Ixᵃᵃᶠ(i, j, k, grid, c) = @inbounds T(0.5) * (c[i, j,   k] - c[i, j, k-1]) 

@inline ∇⋅ᶜᶜᶜ(i, j, k, grid, u, v, w) = (
                                         ∂xᶜᵃᵃ(i, j, k, grid, u) +
                                         ∂yᵃᶜᵃ(i, j, k, grid, v) +
                                         ∂zᵃᵃᶜ(i, j, k, grid, w)
                                        )

#
# Thermodynamics
#

struct DryIdealGas
    R
    cp
    cv
    T0
    ρ0
    s0
end

function diagnose_pressure(gas, ρ, s)
    p0 = gas.ρ0 * gas.R * gas.T0
    return p0 * exp((s - gas.s0 + gas.cp * log(ρ/gas.ρ0)) / gas.cv)
end

#
# Boundary conditions
#

function bcc!(grid, c)
    c[-grid.Hx+1:0,:,:] = c[grid.Nx-grid.Hx:grid.Nx,:,:]
    c[grid.Nx+1:grid.Nx+grid.Hx,:,:] = c[1:grid.Hx,:,:]
    c[:,-grid.Hy+1:0,:] = c[:,grid.Ny-grid.Hy:grid.Ny,:]
    c[:,grid.Ny+1:grid.Ny+grid.Hy,:] = c[:,1:grid.Hy,:]
    c[:,:,-grid.Hz+1:0] = c[:,:,1]
    c[:,:,grid.Nz+1:grid.Nz+grid.Hz] = c[:,:,grid.Nz]
end

function bcu!(grid, c)
    bcc!(grid, c)
end

function bcv!(grid, c)
    bcv!(grid, c)
end

function bcw!(grid, c)
    c[-grid.Hx+1:0,:,:] = c[grid.Nx-grid.Hx:grid.Nx,:,:]
    c[grid.Nx+1:grid.Nx+grid.Hx,:,:] = c[1:grid.Hx,:,:]
    c[:,-grid.Hy+1:0,:] = c[:,grid.Ny-grid.Hy:grid.Ny,:]
    c[:,grid.Ny+1:grid.Ny+grid.Hy,:] = c[:,1:grid.Hy,:]
    c[:,:,-grid.Hz+1:1] = 0.0
    c[:,:,grid.Nz+1:grid.Nz+grid.Hz] = 0.0
end

#
# Model
#

struct LinearModel
    grid
    gas
    ρu
    ρv
    ρw
    ρs
    ρ
    p
    ρsu
    ρsv
    ρsw
end

function time_step!(model, Δt)

    # Calculate pressures and entropy fluxes
    for k = 1:model.grid.Nz
        for j = 1:model.grid.Ny
            for i = 1:model.grid.Nx
                model.p[i,j,k] = diagnose_pressure(gas, model.ρ[i,j,k], model.s[i,j,k])
                model.ρsu[i,j,k] = (
                                    model.ρu[i,j,k] * 
                                    Ixᶠᵃᵃ(i, j, k, model.grid, model.ρs) /
                                    Ixᶠᵃᵃ(i, j, k, model.grid, model.ρ)
                                   )
                model.ρsv[i,j,k] = (
                                    model.ρv[i,j,k] * 
                                    Iyᵃᶠᵃ(i, j, k, model.grid, model.ρs) /
                                    Iyᵃᶠᵃ(i, j, k, model.grid, model.ρ)
                                   )
                model.ρsw[i,j,k] = (
                                    model.ρw[i,j,k] * 
                                    Izᵃᵃᶠ(i, j, k, model.grid, model.ρs) /
                                    Izᵃᵃᶠ(i, j, k, model.grid, model.ρ)
                                   )
            end
        end
    end

    # Apply boundary conditions to pressure and entropy fluxes
    bcc!(model.ρ)
    bcc!(model.ρsu)
    bcc!(model.ρsv)
    bcc!(model.ρsw)

    # Calculate tendencies
    for k = 1:model.grid.Nz
        for j = 1:model.grid.Ny
            for k = 1:model.grid.Nz
                model.ρ[i,j,k] += -Δt * ∇⋅ᶜᶜᶜ(i, j, k, model.grid, 
                                              model.ρu, model.ρv, model.ρw)
                model.ρs[i,j,k] += -Δt * ∇⋅ᶜᶜᶜ(i, j, k, model.grid, 
                                               model.ρsu, model.ρsv, model.ρsw)
                model.ρu[i,j,k] += -∂xᶠᵃᵃ(i, j, k, model.grid, model.p)
                model.ρv[i,j,k] += -∂yᵃᶠᵃ(i, j, k, model.grid, model.p)
                model.ρw[i,j,k] += -∂zᵃᵃᶠ(i, j, k, model.grid, model.p)
            end
        end
    end

    # Apply boundary conditions to prognostic variables
    bcc!(model.ρ)
    bcc!(model.ρs)
    bcu!(model.ρu)
    bcv!(model.ρv)
    bcw!(model.ρw)

end

end module
