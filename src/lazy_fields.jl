import Base: getindex
import Oceananigans.Fields: interior

using Base: @propagate_inbounds

using Oceananigans.Operators: ℑxᶠᵃᵃ, ℑyᵃᶠᵃ, ℑzᵃᵃᶠ
using Oceananigans.Fields

struct LazyPrimitiveField{X, Y, Z, A, G, T, C, D} <: AbstractField{X, Y, Z, A, G, T}
          architecture :: A
                  grid :: G
    conservative_field :: C
               density :: D
end

LazyPrimitiveField(LX, LY, LZ, arch, grid, ρϕ, ρ) =
    LazyPrimitiveField{LX, LY, LZ, typeof(arch), typeof(grid), eltype(grid), typeof(ρϕ), typeof(ρ)}(arch, grid, ρϕ, ρ)

@inline @propagate_inbounds getindex(f::LazyPrimitiveField{Center, Center, Center}, I...) =
    @inbounds f.conservative_field[I...] / f.density[I...]

@inline @propagate_inbounds getindex(f::LazyPrimitiveField{Face, Center, Center}, I...) =
    @inbounds f.conservative_field[I...] / ℑxᶠᵃᵃ(I..., f.grid, f.density)

@inline @propagate_inbounds getindex(f::LazyPrimitiveField{Center, Face, Center}, I...) =
    @inbounds f.conservative_field[I...] / ℑyᵃᶠᵃ(I..., f.grid, f.density)

@inline @propagate_inbounds getindex(f::LazyPrimitiveField{Center, Center, Face}, I...) =
    @inbounds f.conservative_field[I...] / ℑzᵃᵃᶠ(I..., f.grid, f.density)

LazyVelocityFields(arch, grid, ρ, ρũ) =
    (u = LazyPrimitiveField(Face, Center, Center, arch, grid, ρũ.ρu, ρ),
     v = LazyPrimitiveField(Center, Face, Center, arch, grid, ρũ.ρv, ρ),
     w = LazyPrimitiveField(Center, Center, Face, arch, grid, ρũ.ρw, ρ))

function LazyTracerFields(arch, grid, ρ, ρc̃)
    c_names = [filter(c -> c != 'ρ', string(c)) for c in keys(ρc̃)]
    c_names = filter(s -> s != "", c_names) .|> Symbol |> Tuple  # Don't include the ρ tracer.

    c_fields = Tuple(
        LazyPrimitiveField(Center, Center, Center, arch, grid, getproperty(ρc̃, Symbol(:ρ, c)), ρ)
        for c in c_names
    )

    return NamedTuple{c_names}(c_fields)
end

function interior(f::LazyPrimitiveField)
    data = zeros(size(f.conservative_field))
    for k in 1:f.grid.Nz, j in 1:f.grid.Ny, i in 1:f.grid.Nx
        data[i, j, k] = f[i, j, k]
    end
    return data
end
