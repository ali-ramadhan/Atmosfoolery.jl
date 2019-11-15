function implicit_w(i, j, grid; T = Float64)

	# Create arrays for holding diagonals of tridiagonal system
	ldiag = Array{T}(undef, grid.Nz-2)
	cdiag = Array{T}(undef, grid.Nz-1)
	udiag = Array{T}(undef, grid.Nz-2)

	# Fill first row (2 non-zero elements)
	k = 2;
	qₖ = ρ⁰[i,j,k] + δτ * (
		-δx_caa(i, j, k, grid, ρu¹) / Δx(i, j, k, grid)
		- δy_aca(i, j, k, grid, ρv¹) / Δy(i, j, k, grid)
		- β⁰ * δz_aac(i, j, k grid, ρw⁰) / ΔzC(i, j, k, grid)
		+ R
		    
		end
	)

	# Fill middle rows (3 non-zero elements)
	for k = 3:grid.Nz-2

	end

	# Fill last row (2 non-zero elements)
	k = grid.Nz-1

end