# Atmosfoolery.jl

[travis-img]: https://travis-ci.com/thabbott/JULES.jl.svg?branch=master
[travis-url]: https://travis-ci.com/thabbott/JULES.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/ni5ifxwqjkbcofsk/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/thabbott/jules-jl

[codecov-img]: https://codecov.io/gh/thabbott/JULES.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/thabbott/JULES.jl

[![travis][travis-img]][travis-url] [![appveyor][appveyor-img]][appveyor-url] [![codecov][codecov-img]][codecov-url]

 Compressible non-hydrostatic model built on top of Oceananigans.jl so it runs on CPUs and GPUs.

# Description

The compressible model implements the conservative Scheme by Satoh (2003), suitable for compressible non-hydrostatic models with moist processes, and is valid in the limit of a condensable gas/atmosphere with multiple moist species. Two choices of prognostic thermodynamic variables are available, internal energy and entropy, although in the Oceananigans spirit adding new thermodynamic variables is pretty easy!

@thabbott implemented the Satoh (2003) equation set. I initially implemented the Klemp et al. (2007) equation set which is only valid in the limit of a dry atmosphere with trace amounts of moist species.

An RK3 time-stepper is used that is apparently 3rd-order accurate for linear terms but only 2nd-order accurate for non-linear terms. This is because it's developed to allow for acoustic time stepping between RK stages. This split-explicit time-stepping scheme is described by Wicker & Skamarock (2002), Klemp et al. (2007), and Satoh (2003). It's essentially the same one used by the NCAR WRF model (Skamarock et al., 2019). No acoustic time stepper is implemented yet. Explicit acoustic time stepping could make sense for regular Cartesian grids while a vertically implicit acoustic time stepper would make sense for vertically stretched grids (possible with `Oceananigans.Solver.BatchedTridiagonalSolver`).

Building the compressible model on top of Oceananigans.jl has allowed it to run on GPUs and make use of the same operators, grids, Coriolis terms, forcing function and boundary conditions, diagnostics, output writers, higher-order advection schemes, and user-interface niceties. Turbulent diffusivity closures may take more work to integrate and not all of them can be shared as the stress tensor is not traceless when the fluid is compressible (see https://github.com/CliMA/Oceananigans.jl/issues/654 for more discussion).

# Reasons why we may consider adding a compressible model to the Oceananigans.jl ecosystem

1. I see Oceananigans.jl as a general-purpose package for fluid dynamics even though we mostly apply it to ocean problems. With both incompressible and compressible models, Oceananigans.jl would appeal to a larger audience and may be used to investigate a greater range of problems.
2. One potential use of the `CompressibleModel` is to simulate a compressible ocean (with pressure as a prognostic variable) in which sound waves artificially slowed down for practical purposes. There were some discussions around this idea and @johncmarshall54 might still be interested.
3. With PR https://github.com/CliMA/Oceananigans.jl/pull/590, Oceananigans.jl will support distributed parallelism via MPI. While incompressible models (and anelastic models) don't scale that well to many nodes due to the need to solve an elliptic Poission equation globally across all ranks, a compressible model is completely local and might easily scale well on many GPUs. So even if the scalability of incompressible models on many GPUs is dissapointing, we may get a very scalable compressible model almost for free. Actually, the efficiency of the Oceananigans MPI algorithm might be best tested using a compressible model.
4. Since distributed FFTs aren't generally available on GPUs (CuFFT only goes up to 16 GPUs), CUDA-aware MPI for incompressible models might take some time and effort to support once PR https://github.com/CliMA/Oceananigans.jl/pull/590 is merged. However, CUDA-aware MPI should work out of the box for compressible models as there are no FFTs to worry about.
5. Due to the need for a fast pressure solver for incompressible models, we are not considering more general grids beyond the vertically stretched Cartesian grid. The compressible model does not have this limitation and can easily make use of a more general Cartesian grid (stretching in all dimensions).
6. The incompressible model is limited to a certain number of topologies, particularly on the GPU, due to the pressure solver. A compressible model would work with all possible topologies out of the box.
7. Since `CompressibleModel` and `IncompressibleModel` share so much common infrastructure they also share a common user interface by construction which makes it easy to switch between the two. I think this is a valueble feature. Most existing packages do not have compressible/incompressible or ocean/atmosphere capabilities within the same package.
8. Under a shared package and user interface, Oceananigans.jl will allow users to easily switch between simulating compressible and incompressible fluids and might also allow for _fast and friendly_ coupled large-eddy simulation (although the amount of work needed to reach this would be non-trivial).

# Mono-repo vs. multiple packages

I think merging this PR puts the Oceananigans.jl repo in danger of becoming a mono-repo so we should be careful.

One big reason why we haven't kept the compressible model in a separate repo is because we just don't have a good name for it yet.

A potential pathway to multiple packages would be to split out the Oceananigans.jl package into four packages: OceananigansBase.jl, OceananigansIncompressible.jl, OceananigansCompressible.jl, and Oceananigans.jl. I'm not sure which modules would go where but the idea is that users will only have to keep interfacing with Oceananigans.jl.

An added advantage of keeping everything under some Oceananigans.jl umbrella is that the name is getting more well-known (and we have a JOSS paper) so I don't think it makes sense to start a second package with a new name that nobody knows (unless it's a good name!).

That said, I would not be opposed to a mono-repo with a well-defined scope. I actually think that this is a better approach. For example, Oceananigans.jl could provide the `CompressibleModel` and `IncompressibleModel` but ocean-specific modules could live in separate packages. We did this with SeawaterPolynomials.jl and could probably do it with other modules to further limit scope if we decide to pursue this approach. So this is still a multiple packages approach but for ancillary features.

# Tests and validation experiments

We spent quite some time ensuring the `CompressibleModel` can simulate some known atmospheric test cases. Following recent trials and tribulations I also decided to add some simple 1D tests. Here I list the tests but will post a followup comment for each test with a figure or animation. Hopefully together these tests act as a starting point to start believing that the `CompressibleModel` indeed does work as expected.

1. Periodic advection of a square waveform
2. Inviscid Burgers equation developing a shock
3. Shock tube problem (Sod, 1978)
4. Rising thermal bubble with entropy and energy (Wicker & Skamarock, 1998)
5. Rising thermal bubble with 3 different gas species (entropy and energy)
6. Density current (Straka et al., 1993)
7. Dry convection

See comments below for movies and eyeball norms.

The four dry rising thermal bubble simulations are used for regression testing.

# GPU performance benchmarks

Preliminary benchmarks show a 75~80x speedup for large models when comparing a single CPU core to a single Titan V GPU on Tartarus. Not as good as the incompressible model as some of the functions that diagnose temperature and pressure need some optimizing, especially in the case of multiple gases.

```
Compressible model benchmarks                                                                                                                                                                                     
┌──────┬──────┬───────────┬───────────┬────────────┬────────────┬────────────┬────────────┬────────────┬────────┐
│ Arch │ Size │     Gases │ ThermoVar │        min │     median │       mean │        max │     memory │ allocs │
├──────┼──────┼───────────┼───────────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────┤
│  CPU │  32³ │  DryEarth │    Energy │  37.682 ms │  38.012 ms │  37.982 ms │  38.199 ms │ 646.33 KiB │   4499 │
│  CPU │  32³ │  DryEarth │   Entropy │  32.325 ms │  32.920 ms │  32.928 ms │  33.628 ms │ 646.33 KiB │   4499 │
│  CPU │  32³ │ DryEarth3 │    Energy │  52.473 ms │  52.815 ms │  52.896 ms │  53.413 ms │ 816.50 KiB │   5635 │
│  CPU │  32³ │ DryEarth3 │   Entropy │  69.928 ms │  70.388 ms │  70.402 ms │  71.224 ms │ 816.50 KiB │   5635 │
│  CPU │ 192³ │  DryEarth │    Energy │    7.438 s │    7.438 s │    7.438 s │    7.438 s │ 646.33 KiB │   4499 │
│  CPU │ 192³ │  DryEarth │   Entropy │    6.501 s │    6.501 s │    6.501 s │    6.501 s │ 646.58 KiB │   4512 │
│  CPU │ 192³ │ DryEarth3 │    Energy │   11.265 s │   11.265 s │   11.265 s │   11.265 s │ 816.75 KiB │   5648 │
│  CPU │ 192³ │ DryEarth3 │   Entropy │   15.038 s │   15.038 s │   15.038 s │   15.038 s │ 816.50 KiB │   5635 │
│  GPU │  32³ │  DryEarth │    Energy │   8.328 ms │   8.513 ms │   8.608 ms │   9.676 ms │   2.00 MiB │  30129 │
│  GPU │  32³ │  DryEarth │   Entropy │   7.863 ms │   8.500 ms │   8.529 ms │   9.515 ms │   2.00 MiB │  30129 │
│  GPU │  32³ │ DryEarth3 │    Energy │  10.133 ms │  10.754 ms │  10.751 ms │  11.281 ms │   2.54 MiB │  37805 │
│  GPU │  32³ │ DryEarth3 │   Entropy │   9.992 ms │  10.572 ms │  10.542 ms │  10.836 ms │   2.54 MiB │  37839 │
│  GPU │ 192³ │  DryEarth │    Energy │ 101.341 ms │ 101.612 ms │ 101.589 ms │ 101.709 ms │   2.01 MiB │  30270 │
│  GPU │ 192³ │  DryEarth │   Entropy │  86.051 ms │  86.195 ms │  86.226 ms │  86.710 ms │   2.01 MiB │  30270 │
│  GPU │ 192³ │ DryEarth3 │    Energy │ 139.732 ms │ 140.009 ms │ 139.957 ms │ 140.079 ms │   2.54 MiB │  37983 │
│  GPU │ 192³ │ DryEarth3 │   Entropy │ 375.725 ms │ 376.142 ms │ 376.123 ms │ 376.399 ms │   2.54 MiB │  37983 │
└──────┴──────┴───────────┴───────────┴────────────┴────────────┴────────────┴────────────┴────────────┴────────┘
```

```
Compressible model speedups                                                                                                                                                                                       
┌──────┬───────────┬───────────┬─────────┐
│ Size │     Gases │ ThermoVar │ speedup │
├──────┼───────────┼───────────┼─────────┤
│  32³ │  DryEarth │    Energy │  4.465x │
│  32³ │  DryEarth │   Entropy │  3.873x │
│  32³ │ DryEarth3 │    Energy │  4.911x │
│  32³ │ DryEarth3 │   Entropy │  6.658x │
│ 192³ │  DryEarth │    Energy │ 73.203x │
│ 192³ │  DryEarth │   Entropy │ 75.421x │
│ 192³ │ DryEarth3 │    Energy │ 80.457x │
│ 192³ │ DryEarth3 │   Entropy │ 39.981x │
└──────┴───────────┴───────────┴─────────┘
```

# TODO

Right now everything lives in a `compressible` directory to keep it separate.

There are many improvements that could be made to the `CompressibleModel`, starting first with how reference states are initialized and how initial conditions are set (see verification scripts), but this is a list of TODO items in case we decide to merge this PR.

- [ ] Add some documentation, especially for the numerical methods. I should probably LaTeX some notes first [@thabbot has some but might be for the Klemp et al. (2007) equations?].
- [ ] Transfer issues from JULES.jl to Oceananigans.jl to preserve useful discussions and action items.
- [ ] Merge modules, tests, and verification experiments from `compressible` directory.

# References

* Jahn et al. (2015): https://doi.org/10.5194/gmd-8-317-2015
* Klemp et al. (2007): https://doi.org/10.1175/MWR3440.1
* Satoh (2003): https://doi.org/10.1175/1520-0493(2003)131%3C1033:CSFACN%3E2.0.CO;2
* Skamarock et al. (2019): https://opensky.ucar.edu/islandora/object/opensky%3A2898
* Sod (1978): https://doi.org/10.1016/0021-9991(78)90023-2
* Straka et al. (1993): https://doi.org/10.1002/fld.1650170103
* Wicker and Skamarock (1998): https://doi.org/10.1175/1520-0493(1998)126%3C1992:ATSSFT%3E2.0.CO;2

# 1. Periodic advection of a square waveform

Using higher-order advection schemes from Oceananigans.jl to advect momentum and tracers.

![periodic_advection_N64_CFL0 10_CenteredFourthOrder_U-1](https://user-images.githubusercontent.com/20099589/96395228-15964880-1192-11eb-887a-a656049c5d57.gif)

Older figure with WENO-5 (too lazy to re-generate):

![95659750-3b1bb600-0af1-11eb-99f4-fcbe44541bc4](https://user-images.githubusercontent.com/20099589/96395215-116a2b00-1192-11eb-8024-ad0a2a7dd81e.gif)

# 2. Inviscid Burgers equation developing a shock

Was expecting only one shock to form (the one on the left) but probably due to bad setup there's another one on the right that blows up. So probably not the best validation experiment lol.

![burgers_equation](https://user-images.githubusercontent.com/20099589/96395301-4c6c5e80-1192-11eb-9c51-395b67926f74.gif)

# 3. Shock tube

Was too lazy to code up the analytical solution. See Wikipedia for what it should look like: https://en.wikipedia.org/wiki/Sod_shock_tube

With `CenteredSecondOrder`:

![shock_tube_CenteredSecondOrder](https://user-images.githubusercontent.com/20099589/96395425-9d7c5280-1192-11eb-98e6-7434cf011788.gif)

With `WENO5` [little wiggles might be because I'm dividing by ρ instead of WENO5(ρ) or something?]:

![shock_tube](https://user-images.githubusercontent.com/20099589/96395430-a10fd980-1192-11eb-88a1-949e849023a1.gif)

# 4. Inviscid dry rising thermal bubble

WENO-5 doing pretty well with ν = κ = 0, some nice Kelvin-Helmholtz instabilities.

See YouTube: https://www.youtube.com/watch?v=_CgyNHLlpvA

![image](https://user-images.githubusercontent.com/20099589/96395700-504cb080-1193-11eb-8198-ea13031c772a.png)

# 5. Viscous rising thermal bubble with 3 different gas species

A simulation of a warm thermal bubble rising through a viscous compressible fluid made up of three gases with different density profiles. Credit to @thabbott for getting this to work!

See YouTube: https://www.youtube.com/watch?v=mpdRlLd1PTM

![image](https://user-images.githubusercontent.com/20099589/96395848-b1748400-1193-11eb-8e67-32e06ba3aabf.png)

# 6. Nonlinear density current

Basically the rising thermal bubble in reverse.

See YouTube: https://www.youtube.com/watch?v=guk77zjHR6k

![image](https://user-images.githubusercontent.com/20099589/96395940-f26c9880-1193-11eb-8bf3-b7b087f712b5.png)

# 7. Dry convection

Just playing around to see if time-stepping a larger 3D model on a GPU would uncover any issues (it didn't :tada:). And it was an excuse to use Makie.jl.

See YouTube: https://www.youtube.com/watch?v=O0TTfkgNRCI

![image](https://user-images.githubusercontent.com/20099589/96396098-568f5c80-1194-11eb-8c50-f51b6fd4d0f5.png)
