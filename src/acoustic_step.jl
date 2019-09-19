include("acoustic_operators.jl")

function pressure_perturbation(ρ⁰, ρs⁰, pt, ρt, cv, cp)
	return pt ./ (cv .* ρt) .* (ρs⁰ .+ cp .* ρ⁰)
end

function implicit_w(ρu⁰, ρv⁰, ρw⁰, ρ⁰, ρs⁰, grid, δt, Rw, Rρ, Rs, ρt, st, pt, α)

	# Set up storage
	w = WFaceField(grid)
	rhs = Array{eltype(w)}(undef, grid.nz-1)
	udiag = Array{eltype(w)}(undef, grid.nz-2)
	diag = Array{eltype(w)}(undef, grid.nz-1)
	ldiag = Array{eltype(w)}(undef, grid.nz-2)

	# Calculate time differencing coefficients
	β⁰ = (1.0 - α) / 2.0
	β¹ = (1.0 + α) / 2.0

	# Create coefficient matrix
	# TODO: implement
	udiag .= (

	)
	diag .= (
		
	)
	ldiag .= (
		
	)
	coef = Tridiagonal(ldiag, diag, udiag)

	# Loop over columns
	for ii = 1:grid.nx
		for jj = 1:grid.ny

			# Calculate RHS 
			# TODO: implement

			# Solve
			w[2:grid.nz] .= coef \ rhs

		end
	end

end

"""
	acoustic_step(ρu⁰, ρv⁰, ρw⁰, ρ⁰, ρs⁰, p⁰, δt, eos, Ru, Rv, Rw, Rρ, Rs, s)

Integrate perturbation fields for momentum, mass, and entropy over a single acoustic time step using the algorithm described in Klemp et. al. (2007).
"""
function acoustic_step(
	ρu⁰, ρv⁰, ρw⁰, ρ⁰, ρs⁰, # Perturbation fields
	grid, δt, α				# Information about grid and time step
	Ru, Rv, Rw, Rρ, Rs, 	# Large time step tendencies
	ρt, st, pt, 			# Fields at start of large time step
	g, cv, cp,				# Physical constants
	)

	# Calculate perturbation pressure field
	p⁰ = eos(ρ⁰, ρs⁰, pt, ρt, cv, cp)
	# Update ghost points?

	# Advance horizontal momentum
	ρu¹ = ρu⁰ .+ δt .* (-∂x_c2f(p⁰) .+ Ru)
	ρv¹ = ρv⁰ .+ δt .* (-∂y_c2f(p⁰) .+ Rv)
	# Update ghost points?

	# Solve implicit equation for vertical momentum
	ρw¹ = implicit_w(ρu⁰, ρv⁰, ρw⁰, ρ⁰, ρs⁰, grid, δt, Rw, Rρ, Rs, ρt, st, pt, α)


	# Advance mass and entropy
	ρw¹² = avgt(ρw⁰, ρw¹)
	ρ¹ = ρ⁰ .+ δt .* (-divg(ρu¹, ρv¹, ρw¹²) + Rρ)
	su = avgx_c2f(st)
	sv = avgy_c2f(st)
	sw = avgz_c2f(st)
	ρs¹ = ρs⁰ .+ δt .* (-divg(ρu¹, ρv¹, ρw¹²; su = su, sv = sv, sw = sw) + Rs)

	return (ρu¹, ρv¹, ρw¹, ρ¹, ρs¹)

end