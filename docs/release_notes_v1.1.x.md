---
title: Version 1.1 Update Notes
parent: Update Notes
nav_order: 1
---

# Version 1.1 Update Notes

The v1.1 release line covers numerical failure handling, vortex-example
stability, macOS portability, and tutorial validation, preserving the existing
input format and solver interfaces wherever possible.

## Version 1.1

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
- Homebrew GNU Fortran is supported without the Linux-only `-mcmodel=medium`
  flag, and external libraries can be built from any working directory with
  explicit compiler choices, as in
  `FC=gfortran CC=cc BUILD_JOBS=2 ./external/CMake_build.sh`.
- Documented and notebook MPI commands use portable `mpirun -n`, and the
  validation tool adds OpenMPI-only flags at runtime only when they are
  required, so compatible MPICH installations can also run the tutorials.
- Added the tutorial validation gate, which checks six notebooks, twelve
  documentation pages, and 13 documented executable paths on the requested MPI
  rank count, runs matched low-viscosity SVV-off and SVV-on vortex cases and
  requires their final fields to differ, and requires representative
  zero-coefficient hyperviscosity, unsupported-power, and invalid-SVV inputs to
  fail nonzero. Notebook cells are executed as well when `nbformat`, `nbclient`,
  and a Python kernel are installed. CI repeats the gate with GNU Fortran bounds
  checking, backtraces, and floating-point traps enabled.
- Known limitation: exact 2/3 de-aliasing applies to the periodic azimuthal and
  axial directions only. The mapped radial quadratic product is not
  over-integrated; SVV can damp its high-mode tail but cannot undo radial energy
  already aliased into retained modes.
- Known limitation: the saved vortex sample field is vorticity magnitude, not a
  kinetic-energy spectrum. A velocity-spectrum diagnostic and resolution study
  are required before assessing a `-5/3` inertial-range slope.
- Known limitation: with periodic 2/3 masking, `svv_cutoff = 0.75` primarily
  targets the mapped radial tail. Use a cutoff below approximately `2/3` if the
  highest retained periodic bands should also participate, and retune the target
  and strength for the chosen Reynolds number and resolution.
- Known issue, fixed in v1.1.1: the documented gate command
  `python3 tools/validate_tutorials.py --ranks 2` always failed on this version,
  because `wave_propagation_1d` stops unless given exactly one process and the
  gate passed its global rank count to every application. On v1.1 itself, use
  `--ranks 1`.

## Version 1.1.1

- Fixed the tutorial validation gate, which could not pass above one MPI rank.
  `wave_propagation_1d` solves a purely radial problem, so the rank owning the
  radial direction must own all of it, and the program stops deliberately when
  launched on more than one process. The gate passed its global `--ranks` value
  to every application, so that case aborted at any rank count above one. The
  rank count is now overridden for that application alone, and the twelve
  remaining tutorials keep their parallel coverage.
- Corrected `CITATION.cff`, which recorded the v1.0 Zenodo DOI together with a
  version string and release date matching no tag.
- This maintenance patch did not change the numerical solver, input format, or
  documented tutorial commands; the sources under `src/` are identical to v1.1.

## Version 1.1.2

- Combined the 1.1 patch update notes into this single page, matching the
  `v1.0.x` convention of one page per major-minor line, and aligned its
  structure and wording with the v1.0 page.
- Updated the validation gate's documentation inventory and the Update Notes
  index to match the merged page.
- This maintenance patch did not change the numerical solver, input format, or
  documented tutorial commands; the sources under `src/` are identical to v1.1.

## Version 1.1.3

- `CITATION.cff` now records the concept DOI rather than a version DOI. A
  version DOI is minted from the release it identifies, so no tag can contain
  its own; the concept DOI always resolves to the newest version and is correct
  in every release. The version DOIs for v1.1, v1.1.1 and v1.1.2 are listed
  under `identifiers`.
- This maintenance patch did not change the numerical solver, input format, or
  documented tutorial commands; the sources under `src/` are identical to v1.1.
