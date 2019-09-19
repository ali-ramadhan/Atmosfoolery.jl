"""
	∂x_c2f(s, grid)

Calculate the x derivative of a cell-centered field on cell faces.
"""
function ∂x_c2f(s::CellField, grid::Grid)
	u = UFaceField(grid)
	fu[1:grid.nz,1:grid.nx,1:grid.nx+1] .= (
		s[1:grid.nz,1:grid.ny,1:grid.nx+1] .- 
		s[1:grid.nz,1:grid.ny,0:grid.nx]
	) ./ grid.δx
	return u
end

"""
	∂y_c2f(s, grid)

Calculate the y derivative of a cell-centered field on cell faces.
"""
function ∂y_c2f(s::CellField, grid::Grid)
	v = VFaceField(grid)
	v[1:grid.nz,1:grid.ny+1,1:grid.nx] .= (
		s[1:grid.nz,1:grid.ny+1,1:grid.nx] .- 
		s[1:grid.nz,0:grid.ny,1:grid.nx]
	) ./ grid.δy
	return v
end

"""
	∂z_c2f(s, grid)

Calculate the z derivative of a cell-centered field on cell faces.
"""
function ∂z_c2f(s::CellField, grid::Grid)
	w = WFaceField(grid)
	w[1:grid.nz+1,1:grid.ny,1:grid.nx] .= (
		s[1:grid.nz+1,1:grid.ny,1:grid.nx] .- 
		s[0:grid.nz,1:grid.ny,1:grid.nx]
	) ./ grid.δzf
	return w
end

"""
	∂z_f2c(w, grid)

Calculate the z derivative of a face-centered field on cell centers
"""
function ∂z_f2c(w::WFaceField, grid::Grid)
	s = CellField(grid)
	s[1:grid.nz,1:grid.ny,1:grid.nx] .= (
		w[2:grid.nz+1,1:grid.ny,1:grid.nx] .- 
		w[1:grid.nz,1:grid.ny,1:grid.nx]
	) ./ grid.δzc
	return s
end

"""
	avgx_c2f(s, grid)

Interpolate cell-centered fields onto cell faces in the x direction.
"""
function avgx_c2f(s::CellField, grid::Grid)
	u = UFaceField(grid)
	u[1:grid.nz,1:grid.ny,1:grid.nx+1] .= 0.5 .* (
		s[1:grid.nz,1:grid.ny,1:grid.nx+1] .+ s[1:grid.nz,1:grid.ny,0:grid.nx]
	)
	return u
end

"""
	avgy_c2f(s, grid)

Interpolate cell-centered fields onto cell faces in the y direction.
"""
function avgy_c2f(s::CellField, grid::Grid)
	v = VFaceField(grid)
	v[1:grid.nz,1:grid.ny+1,1:grid.nx] .= 0.5 .* (
		s[1:grid.nz,1:grid.ny+1,1:grid.nx] + s[1:grid.nz,0:grid.ny,1:grid.nx]
	)
	return v
end

"""
	avgz_c2f(s, grid)

Interpolate cell-centered fields onto cell faces in the z direction.
"""

function avgz_c2f(s::CellField, grid::Grid)
	w = WFaceField(grid)
	w[1:grid.nz+1,1:grid.ny,1:grid.nx] .= 0.5 .* (
		s[1:grid.nz+1,1:grid.ny,1:grid.nx] + s[0:grid.nz,1:grid.ny,1:grid.nx]
	)
	return w
end

"""
	divg(ρu, ρv, ρw, su, sv, sw)

Calculate the tracer flux out of a cell produced by cell face mass fluxes ρu, ρv, and ρw carrying tracers with densities su, sv, and sw.
"""
function divg(ρu::UFaceField, ρv::VFaceField, ρw::WFaceField, grid::Grid; 
	su = 1.0, sv = 1.0, sw = 1.0)
	
	# Create field to hold result
	s = CellField(grid)

	# Calculate fluxes
	Fu = ρu .* su
	Fv = ρv .* sv
	Fw = ρw .* sw

	# Calculate net flux into cell in each direction
	divFu = grid.Au .* (
		Fu[1:grid.nz,1:grid.ny,2:grid.nx+1] .- Fu[1:grid.nz,1:grid.ny,1:grid.nx]
	)
	divFv = grid.Av .* (
		Fv[1:grid.nz,2:grid.ny+1,1:grid.nx] .- Fu[1:grid.nz,1:grid.ny,1:grid.nx]
	)
	divFw = grid.Aw .* (
		Fw[2:grid.nz+1,1:grid.ny,1:grid.nx] .- Fw[1:grid.nz,1:grid.ny,1:grid.nx]
	)

	# Normalize by cell volume and return
	s[1:grid.nz,1:grid.ny,1:grid.nx] .= (divFu .+ divFv .+ divFw) ./ grid.V
	return s

end