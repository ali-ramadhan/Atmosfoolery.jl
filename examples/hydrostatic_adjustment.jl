using JULES

# Create grid
grid = Grid(4, 4, 4, 1, 1, 1, 10.0, 10.0, 10.0)

# Create prognostic fields
ρ = CellField(grid)
ρs = CellField(grid)
ρu = FaceFieldX(grid)
ρv = FaceFieldY(grid)
ρw = FaceFieldZ(grid)

# Create diagnostic fields
p = CellField(grid)
ρsu = FaceFieldX(grid)
ρsv = FaceFieldY(grid)
ρsw = FaceFieldZ(grid)

# Define gas constants
T0 = 273.16
p0 = 1e5
R = 287.0
ρ0 = p0 / (R*T0)
gas = DryIdealGas(R, 3.5*R, 2.5*R, T0, ρ0, 0.0)

# Define planetary constants
world = World(9.81)

# Create model
model = LinearModel(grid, gas, world,
                    ρu, ρv, ρw, ρs, ρ,
                    p, ρsu, ρsv, ρsw)

# Initialize model
initialize_isothermal_atmosphere!(model, 300.0, 1e5)

# Try time stepping
for _ = 1:1
    time_step!(model, 1e-6)
end
