---
title: Version 1.1 Update Notes
parent: Update Notes
nav_order: 2
---

# Version 1.1 Update Notes

Version 1.1 improves numerical failure handling, vortex-example stability,
macOS portability, and tutorial validation while preserving the existing input
format and solver interfaces wherever possible.

## Numerical and solver corrections

- The 3D vortex example recomputes its nonlinear right-hand side after startup
  extrapolation and uses a local final-step size without modifying the global
  time step.
- Periodic quadratic products use a 2/3 Fourier mask. An optional bounded
  tail-energy spectral vanishing-viscosity (SVV) filter can damp populated high
  modes; separate feedback states are maintained for the toroidal and poloidal
  potentials.
- The mapped far-field coefficient is removed after each vortex step, and all
  MPI ranks collectively abort if either potential becomes non-finite.
- Laplacian inversion null modes, transformed-input handling, FFT chopped-region
  checks, and dense-by-band storage bounds are corrected.
- Invalid chopping sizes, unsupported hyperviscosity powers, enabled
  hyperviscosity with a non-positive coefficient, malformed SVV controls, and
  invalid output intervals now fail during input parsing.

## macOS and MPI portability

Homebrew GNU Fortran is supported without the Linux-only `-mcmodel=medium`
flag. External libraries can be built from any working directory with explicit
compiler choices:

```bash
FC=gfortran CC=cc BUILD_JOBS=2 ./external/CMake_build.sh
make
```

Documented and notebook MPI commands use portable `mpirun -n`; the validation
tool adds OpenMPI-only flags at runtime when they are actually required. This
also allows compatible MPICH installations to run the tutorials.

## Validation gate

The standard gate validates six notebooks, twelve documentation pages, and 13
documented executable paths on the requested MPI rank count. It additionally
runs matched low-viscosity SVV-off and SVV-on vortex cases and requires their
final fields to differ. Representative zero-coefficient hyperviscosity,
unsupported-power, and invalid-SVV inputs must also fail nonzero.

```bash
python3 tools/validate_tutorials.py --ranks 2
```

{: .warning }
> In v1.1 this command always failed. `wave_propagation_1d` is a purely radial
> program that stops unless it is given exactly one process, and the gate passed
> its global `--ranks` value to every application, so the case aborted at any
> rank count above one. Fixed in
> [v1.1.1]({{ site.baseurl }}/release_notes_v1.1.1.html), which overrides the
> rank count for that application only. On v1.1, use `--ranks 1`.

Notebook cells can also be executed when `nbformat`, `nbclient`, and a Python
kernel are installed:

```bash
python3 tools/validate_tutorials.py --execute-notebooks --ranks 1
```

CI repeats the two-rank gate with GNU Fortran bounds checking, backtraces, and
floating-point traps enabled.

## Numerical limitations and tuning

- The current exact 2/3 de-aliasing applies to the periodic azimuthal and axial
  directions. The mapped radial quadratic product is not over-integrated; SVV
  can damp its high-mode tail but cannot undo radial energy already aliased into
  retained modes.
- The saved vortex sample field is vorticity magnitude, not a kinetic-energy
  spectrum. A velocity-spectrum diagnostic and resolution study are required
  before assessing a `-5/3` inertial-range slope.
- With periodic 2/3 masking, `svv_cutoff = 0.75` primarily targets the mapped
  radial tail. Use a cutoff below approximately `2/3` if the highest retained
  periodic bands should also participate, and retune the target and strength
  for the chosen Reynolds number and resolution.
