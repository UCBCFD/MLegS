---
title: Version 1.0 Update Notes
parent: Update Notes
nav_order: 2
---

# Version 1.0 Update Notes

The v1.0 release line covers the initial open-source MLegS package and its first
metadata and publication maintenance updates.

## Version 1.0

- Established the modernized and parallelized MLegS codebase for mapped
  Legendre spectral simulations in radially unbounded domains.
- Introduced the current Modern Fortran modules, MPI-enabled distributed
  scalar operations, spectral transforms, PDE time integrators, vector-field
  utilities, and executable sample programs.
- Added the tutorial notebooks and Jekyll documentation for prerequisites,
  spectral transformations and operations, time integration, vector fields,
  and representative scalar and vortical-flow problems.
- Added the CMake-based external-library build flow for FFTE, FM, and LAPACK,
  together with Makefile targets for modules and applications.
- Moved the GitHub Pages deployment to the `main` branch so the documentation
  follows the open-source release line.

## Version 1.0.1

- Added and standardized `CITATION.cff` with the software title, authors,
  affiliations, ORCID, BSD-2-Clause license, release metadata, and Zenodo DOI
  `10.5281/zenodo.14976471`.
- Added the README DOI badge and project image.
- This maintenance patch did not change the numerical solver, input format,
  or documented tutorial commands.

## Version 1.0.2

- Updated the post-open-source publication entry for transient wake-vortex
  growth from its arXiv record to the published *Journal of Fluid Mechanics*
  article, volume 1014, A16, with DOI `10.1017/jfm.2025.253`.
- This maintenance patch did not change the numerical solver, input format,
  or tutorial execution paths.
