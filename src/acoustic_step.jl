include("acoustic_operators.jl")

"""
	acoustic_step(ρu⁰, ρv⁰, ρw⁰, ρ⁰, ρs⁰, p⁰, δt, eos, Ru, Rv, Rw, Rρ, Rs, s)

Integrate perturbation fields for momentum, mass, and entropy over a single acoustic time step using the algorithm described in Klemp et. al. (2007).
"""
function acoustic_step(ρu⁰, ρv⁰, ρw⁰, ρ⁰, ρs⁰, δt, eos, Ru, Rv, Rw, Rρ, Rs, st)

	# Advance horizontal momenta
	p⁰ = eos(ρ⁰, ρs⁰)
	ρu¹ = ρu⁰ .+ δt .* (-∂x_c2f(p⁰) .+ Ru)
	ρv¹ = ρv⁰ .+ δt .* (-∂y_c2f(p⁰) .+ Rv)

	# Solve implicit equation for vertical momentum


	# Advance mass and entropy
	ρw¹² = avgt(ρw⁰, ρw¹)
	ρ¹ = ρ⁰ .+ δt .* (-divg(ρu¹, ρv¹, ρw¹²) + Rρ)
	su = avgx_c2f(st)
	sv = avgy_c2f(st)
	sw = avgz_c2f(st)
	ρs¹ = ρs⁰ .+ δt .* (-divg(ρu¹, ρv¹, ρw¹²; su = su, sv = sv, sw = sw) + Rs)

	return (ρu¹, ρv¹, ρw¹, ρ¹, ρs¹)

end