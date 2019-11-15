struct Grid 
	nx::Int64
	ny::Int64
	nz::Int64
	δx::Float64
	δy::Float64
	δzf::OffsetVector{Float64}
	δzc::OffsetVector{Float64}
	Au::OffsetVector{Float64}
	Av::OffsetVector{Float64}
	Aw::Float64
	V::OffsetVector{Float64}
end

function isotropic_grid(; n = 4, δ = 10.0)
	nx = n
	ny = n
	nz = n
	δx = δ
	δy = δ
	δzf = OffsetVector{Float64}(δ, 1:nz+1)
	δzc = OffsetVector{Float64}(δ, 1:nz)
	Au = δy .* δzc
	Av = δx .* δzc
	Aw = δx .* δy
	V = δx .* δy .* δzc
	return Grid(nx, ny, nz, δx, δy, δzf, δzc, Au, Av, Aw, V)
end

const CellField = OffsetArray{Float64}
const UFaceField = OffsetArray{Float64}
const VFaceField = OffsetArray{Float64}
const WFaceField = OffsetArray{Float64}

function cellfield(grid::Grid)
	return CellField(undef, 0:grid.nz+1, 0:grid.ny+1, 0:grid.nx+1)
end

function ufacefield(grid::Grid)
	return UFaceField(undef, 1:grid.nz, 1:grid.ny, 1:grid.nx+1)
end

function vfacefield(grid::Grid)
	return VFaceField(undef, 1:grid.nz, 1:grid.ny+1, 1:grid.nx)
end

function wfacefield(grid::Grid)
	return WFaceField(undef, 1:grid.nz+1, 1:grid.ny, 1:grid.nx)
end