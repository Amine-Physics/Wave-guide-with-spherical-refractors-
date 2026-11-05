# Waveguide with Circular Refractors: a 2D Finite-Element Scattering Solver

A 2D frequency-domain solver for a scalar wave in a rectangular waveguide containing one or several circular dielectric refractors. It is built with [deal.II](https://dealii.org) and meshes with the [Gmsh](https://gmsh.info) C++ API. The long-term goal is to study bound states in the continuum (BICs); the current version is a validated scattering solver with mode-resolved post-processing.

## Physics

The unknown `u` is the out-of-plane magnetic field (TM case) and satisfies

    div( (1/eps_r) grad u ) + k0^2 u = 0,        k0 = omega / c

inside the guide of width `W = y_max - y_min`, where `eps_r` is 1 in the background and constant inside each refractor.

- **Walls** (boundary id 3): `u = 0`.
- **Ports** (boundary ids 1 = input, 2 = output): a first-order Robin condition built with the propagation constant of the incident mode `n`,

      beta_n = sqrt( k0^2 - (n*pi/W)^2 ).

  It absorbs mode `n` exactly. At the input port a source term injects a unit-normalized incident wave. In the code's convention the forward wave varies as `exp(-i*beta*x)`.
- **Single-mode window:** only mode 1 propagates when `pi/W < k0 < 2*pi/W`. Outside it, several modes propagate and the single-beta port is only exact for the incident mode (see Limitations).

The complex field is stored as two real components (`Re`, `Im`) with `FESystem<2>(FE_Q<2>(1)^2)`, and the linear system is solved with UMFPACK.

## Mesh

`build_mesh` (Gmsh C++ API, OpenCASCADE kernel) builds the rectangle and the disks, makes the disk interfaces conforming with a boolean fragment, and writes a quadrilateral mesh (`msh` v2.2). Labels carried into deal.II:

| Label | Meaning |
|---|---|
| material id 1 | background |
| material id 2 + i | refractor `i` (same order as in the parameter file) |
| boundary id 1 / 2 | input port / output port |
| boundary id 3 | walls |

After reading the mesh, each refractor gets its own `SphericalManifold` (manifold id = refractor index), set on the interface faces and the adjacent cells, before global refinement. Keep gaps between refractors (and to the walls) several cells wide, otherwise a coarse cell can touch two interfaces.

## Post-processing

- **Modal amplitude** `a_m(x)`: projection of the field on the mode `sin(m*pi*(y - y_min)/W)` along the line `x = const`, using Gauss quadrature and `FEFieldFunction`.
- **Forward/backward fit:** at two planes on each side of the refractors, `a_m(x) = A*exp(-i*beta_m*x) + B*exp(+i*beta_m*x)` is solved for `A` and `B`.
- **Scattering results:** `r = B_left/A_left`, `t = A_right/A_left`, `R = |r|^2`, `T = |t|^2`, balance `R + T`.
- **Power spectrum:** for every propagating mode `m`, the powers `beta_m*|B_left|^2`, `beta_m*|A_right|^2` and the port-reflection terms, normalized by the incident power. The sum over modes is the generalized energy balance.
- **Intensity:** transverse intensity `integral |u|^2 dy` along the guide.

## Build

Requirements: CMake, a C++17 compiler, deal.II (>= 9.3, with MPI), and the Gmsh C++ library and headers (e.g. `libgmsh-dev`).

```
mkdir builds && cd builds
cmake ..
make
```

If CMake cannot find Gmsh, pass its location, e.g. `cmake -DCMAKE_PREFIX_PATH=/path/to/gmsh ..`. Run from `build/`, with the parameter file in the same directory:

```
./main        # or: make run
```

## Parameters

The parameter file is read once into a plain `Parameters` struct (`read_parameters`). Subsections and keys:

```
subsection Mesh & Geometry
  set Number Of Refinements = 2
  set Wave Guide Corners    = -1.0, 0.0, 1.0, 0.5   # x_min, y_min, x_max, y_max
  set Max Cell Size         = 0.05
end

subsection Refractors            # lists of equal length, one entry per refractor
  set Centers x      = 0.0, 0.6
  set Centers y      = 0.25, 0.25
  set Radii          = 0.15, 0.10
  set Permittivities = 1.5, 2.0
end

subsection Physics
  set Omega       = 3e9      # angular frequency in rad/s (f = omega / 2 pi)
  set C           = 3e8
  set Mode Number = 1
end

subsection Output
  set Output File Name = solution
end
```

Post-processing options (number of quadrature subintervals, sampling points) are declared in the same file; see `parameters.cc`. Refractors must lie inside the guide and must not overlap; both are checked on reading.

## Output

- `<name>_<n>.vtu` and `.h5`: fields `Re`, `IM` and `material_id`, for ParaView.
- `grid.svg`: the refined mesh.
- CSV files: modal amplitude and intensity along the guide (sampled points and a uniform graph), `scattering_results.csv`, and the per-mode power spectrum.

The scripts in `post_processing/` read the CSV files and produce the amplitude, phase, intensity and power-spectrum figures.

## Repository layout

```
CMakeLists.txt
main.cc
includes/          headers (Parameters, mesh builder, solver class, post-processor)
src/               implementations
tests/             standalone mesh test (Gmsh only)
post_processing/   Python plotting scripts
```

## Validation status

- Empty guide: reflection close to 0, |t| close to 1, incident amplitude equal to `sqrt(2/W)`.
- Several refractors in the single-mode window: `R + T` close to 1.
- Several propagating modes: the single-mode balance is not expected to equal 1; the per-mode power spectrum is the meaningful output.

## Limitations

- The Robin port is exact only for the incident mode. When other modes propagate, they partly reflect at the ports, which distorts the energy balance (visible as output-port reflection growing for modes near cutoff).
- Q1 elements need roughly 10 or more cells per local wavelength (wavelength divided by `sqrt(eps_r)` inside refractors). High frequencies need fine meshes; check convergence under refinement.
- Measurement planes must be far enough from the refractors for evanescent modes to decay; check that `r` and `t` do not depend on plane positions.
- Only the TM polarization is implemented. For TE, `eps_r` moves from the gradient term to the mass term.

## Roadmap

- Multimode port condition (Dirichlet-to-Neumann map over the propagating modes).
- Convergence study and a comparison with a published waveguide benchmark.
- Eigenvalue / quasi-normal-mode solves with PML or an equivalent open boundary.
- Quality-factor tracking against a symmetry-breaking parameter, to locate symmetry-protected BICs.
- Optional local mesh refinement near the refractors (through Gmsh size fields).

## Author

Amine. Work in progress, developed as preparation for PhD research on electromagnetic simulation. License: to be decided.