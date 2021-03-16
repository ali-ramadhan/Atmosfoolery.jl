"""
This example sets up a cold bubble perturbation which develops into a non-linear
density current. This numerical test case is described by Straka et al. (1993).
Also see: http://www2.mmm.ucar.edu/projects/srnwp_tests/density/density.html

Straka et al. (1993). "Numerical Solutions of a Nonlinear Density-Current -
    A Benchmark Solution and Comparisons." International Journal for Numerical
    Methods in Fluids 17, pp. 1-22.
"""

using Logging
using Printf
using Statistics
using NCDatasets
using CUDA

using Oceananigans
using Oceananigans.Grids
using Oceananigans.Advection
using Oceananigans.OutputWriters
using Oceananigans.Utils
using JULES

using Oceananigans.Fields: cpudata

Logging.global_logger(OceananigansLogger())

const km = kilometers
const hPa = 100.0

Lx = 51.2km
Lz = 6.4km

Δ = 20.0  # grid spacing [m]

Nx = Int(Lx/Δ)
Ny = 1
Nz = Int(Lz/Δ)

topo = (Periodic, Periodic, Bounded)
domain = (x=(-Lx/2, Lx/2), y=(-Lx/2, Lx/2), z=(0, Lz))
grid = RegularCartesianGrid(topology=topo, size=(Nx, Ny, Nz), halo=(3, 3, 3); domain...)

tvar = Energy()

model = CompressibleModel(
                      grid = grid,
                     gases = DryEarth(),
 	         advection = WENO5(),
    thermodynamic_variable = tvar,
                   closure = IsotropicDiffusivity(ν=75.0, κ=75.0)
)

gas = model.gases.ρ
R, cₚ, cᵥ = gas.R, gas.cₚ, gas.cᵥ
g  = model.gravity
pₛ = 1000hPa
Tₛ = 300.0

# Define an approximately hydrostatic background state
θ₀(x, y, z) = Tₛ
p₀(x, y, z) = pₛ * (1 - g*z / (cₚ*Tₛ))^(cₚ/R)
T₀(x, y, z) = Tₛ * (p₀(x, y, z)/pₛ)^(R/cₚ)
ρ₀(x, y, z) = p₀(x, y, z) / (R*T₀(x, y, z))

# Define both energy and entropy
uᵣ, Tᵣ, ρᵣ, sᵣ = gas.u₀, gas.T₀, gas.ρ₀, gas.s₀  # Reference values
ρe₀(x, y, z) = ρ₀(x, y, z) * (uᵣ + cᵥ * (T₀(x, y, z) - Tᵣ) + g*z)
ρs₀(x, y, z) = ρ₀(x, y, z) * (sᵣ + cᵥ * log(T₀(x, y, z)/Tᵣ) - R * log(ρ₀(x, y, z)/ρᵣ))

# Define the initial density perturbation
θᶜ′ = -15.0
xᶜ, zᶜ = 0km, 2km
xʳ, zʳ = 2km, 2km
L(x, y, z) = sqrt(((x - xᶜ)/xʳ)^2 + ((z - zᶜ)/zʳ)^2)
θ′(x, y, z) = (L(x, y, z) <= 1) * θᶜ′ * (1 + cos(π*L(x, y, z))) / 2
ρ′(x, y, z) = -ρ₀(x, y, z) * θ′(x, y, z) / θ₀(x, y, z)

# Define initial state
ρᵢ(x, y, z) = ρ₀(x, y, z) + ρ′(x, y, z)
pᵢ(x, y, z) = p₀(x, y, z)
Tᵢ(x, y, z) = pᵢ(x, y, z) / (R * ρᵢ(x, y, z))

ρeᵢ(x, y, z) = ρᵢ(x, y, z) * (uᵣ + cᵥ * (Tᵢ(x, y, z) - Tᵣ) + g*z)
ρsᵢ(x, y, z) = ρᵢ(x, y, z) * (sᵣ + cᵥ * log(Tᵢ(x, y, z)/Tᵣ) - R * log(ρᵢ(x, y, z)/ρᵣ))

# Set initial state (which includes the thermal perturbation)
set!(model.tracers.ρ, ρᵢ)
tvar isa Energy  && set!(model.tracers.ρe, ρeᵢ)
tvar isa Entropy && set!(model.tracers.ρs, ρsᵢ)
update_total_density!(model)

function print_progress(simulation)
    model, Δt = simulation.model, simulation.Δt
    tvar = model.thermodynamic_variable
    ρᵢ, ρeᵢ, ρsᵢ = simulation.parameters

    zC = znodes(Cell, model.grid)
    ρ̄ᵢ = mean(ρᵢ.(0, 0, zC))
    ρ̄ = mean(cpudata(model.total_density))

    progress = 100 * model.clock.time / simulation.stop_time
    message = @sprintf("[%05.2f%%] iteration = %d, time = %s, CFL = %.4e, acoustic CFL = %.4e, ρ̄ = %.4e (relΔ = %.4e)",
                       progress, model.clock.iteration, prettytime(model.clock.time), cfl(model, Δt),
                       acoustic_cfl(model, Δt), ρ̄, (ρ̄ - ρ̄ᵢ) / ρ̄)

    if tvar isa Energy
        ρ̄ēᵢ = mean(ρeᵢ.(0, 0, zC))
        ρ̄ē = mean(cpudata(model.tracers.ρe))
        message *= @sprintf(", ρ̄ē = %.4e (relΔ = %.4e)", ρ̄ē, (ρ̄ē - ρ̄ēᵢ)/ρ̄ē)
    elseif tvar isa Entropy
        ρ̄s̄ᵢ = mean(ρsᵢ.(0, 0, zC))
        ρ̄s̄ = mean(cpudata(model.tracers.ρs))
        message *= @sprintf(", ρ̄s̄ = %.4e (relΔ = %.4e)", ρ̄s̄, (ρ̄s̄ - ρ̄s̄ᵢ)/ρ̄s̄)
    end

    @info message

    return nothing
end

simulation = Simulation(model, Δt=0.02, stop_time=1000, iteration_interval=50,
                        progress=print_progress, parameters=(ρᵢ, ρeᵢ, ρsᵢ))

fields = Dict(
    "ρ"  => model.total_density,
    "ρu" => model.momenta.ρu,
    "ρw" => model.momenta.ρw
)

tvar isa Energy  && push!(fields, "ρe" => model.tracers.ρe)
tvar isa Entropy && push!(fields, "ρs" => model.tracers.ρs)

simulation.output_writers[:fields] =
    NetCDFOutputWriter(model, fields, filepath="nonlinear_density_current_$(typeof(tvar)).nc",
                       schedule=TimeInterval(10seconds))


# Save base state to NetCDF.
ds = simulation.output_writers[:fields].dataset
ds_ρ = defVar(ds, "ρ₀", Float32, ("xC", "yC", "zC"))
ds_ρe = defVar(ds, "ρe₀", Float32, ("xC", "yC", "zC"))

x, y, z = nodes((Cell, Cell, Cell), grid, reshape=true)
ds_ρ[:, :, :] = ρ₀.(x, y, z)
ds_ρe[:, :, :] = ρe₀.(x, y, z)

run!(simulation)
