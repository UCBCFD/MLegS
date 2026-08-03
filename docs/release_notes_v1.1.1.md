---
title: Version 1.1.1 Update Notes
parent: Update Notes
nav_order: 1
---

# Version 1.1.1 Update Notes

Version 1.1.1 is a maintenance release. It contains no changes to the solver,
the numerical formulation, or any published result; the library sources under
`src/` are identical to v1.1.

## Tutorial validation gate

- `tools/validate_tutorials.py` no longer fails when run with more than one MPI
  rank. `wave_propagation_1d` solves a purely radial problem, so the rank owning
  the radial direction must own all of it, and the program stops deliberately if
  launched on more than one process. The gate previously passed its global
  `--ranks` value to every application, so the documented command
  `python3 tools/validate_tutorials.py --ranks 2` always failed on that one case.
  The rank count is now overridden per application, and the twelve remaining
  tutorials keep their parallel coverage.

## Citation metadata

- `CITATION.cff` records the correct archived-release DOI. The previous entry
  pointed at the v1.0 Zenodo record and carried a version string and date that
  did not match the tag.
