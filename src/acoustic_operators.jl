struct Grid 
	nx
	ny
	nz
	δx
	δy
	δzf
	δzc
end

function isotropic_grid(; n = 4, δ = 10.0)
	nx = n
	ny = n
	nzc = n
	nzf = n+1
	δx = δ
	δy = δ
	δzf = OffsetArray{Float64}(δ, 1:nz+1)
	δzc = OffsetArray{Float64}(δ, 1:nz)
	return Grid(nx, ny, nz, δx, δy, δzf, δzc)
end

"""
	∂x_c2f(s, grid)

Calculate the x derivative of a cell-centered field on cell faces.
"""
function ∂x_c2f(s, grid::Grid)
	return (s[:,:,1:grid.nx+1] .- s[:,:,0:grid.nx]) ./ grid.δx
end

"""
	∂y_c2f(s, grid)

Calculate the y derivative of a cell-centered field on cell faces.
"""
function ∂y_c2f(s, grid::Grid)
	return (s[:,1:grid.ny+1,:] .- s[:,0:grid.ny,:]) ./ grid.δy
end

"""
	∂z_c2f(s, grid)

Calculate the z derivative of a cell-centered field on cell faces.
"""
function ∂z_c2f(s, grid::Grid)
	return (s[1:grid.nzc+1,:,:] .- s[0:grid.nzc,:,:]) ./ grid.δzf
end

"""
	∂z_f2c(w, grid)

Calculate the z derivative of a face-centered field on cell centers
"""
function ∂z_f2c(w, grid::Grid)
	return (w[2:grid.nzf,:,:] .- w[1:grid.nzf-1,:,:]) ./ grid.δzc
end

"""
	avgx_c2f(s, grid)

Interpolate cell-centered fields onto cell faces in the x direction.
"""
function avgx_c2f(s, grid)

end

"""
	avgy_c2f(s, grid)

Interpolate cell-centered fields onto cell faces in the y direction.
"""
function avgy_c2f(s, grid)

end

"""
	avgz_c2f(s, grid)

Interpolate cell-centered fields onto cell faces in the x direction.
"""

function avgz_c2f(s, grid)

end

"""
	divg(ρu, ρv, ρw, su, sv, sw)

Calculate the tracer flux out of a cell produced by cell face mass fluxes ρu, ρv, and ρw carrying tracers with densities su, sv, and sw.
"""
function divg(ρu, ρv, ρw; su = 1.0, sv = 1.0, sw = 1.0)

end