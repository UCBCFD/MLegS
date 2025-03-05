# **MLegS**: Modernized and Parallelized **M**apped **Leg**endre **S**pectral Method Code
[![DOI](https://zenodo.org/badge/836085630.svg)](https://doi.org/10.5281/zenodo.14976470)

MLegS (**M**apped **Leg**endre **S**pectral Method Code) is a code package based on a modernized and parallelized spectral method for numerical simulations in a radially unbounded domain.

<p align="center">
  <img src="https://github.com/user-attachments/assets/c43771b2-e8df-4bd5-b6a4-877a13f839cb" style="border-radius:5px;"/>
</p>

Based on the numerical algorithm proposed by T. Matsushima and P. S. Marcus (1997)[^1], MLegS incorporates scalable multiprocessing interfaces for high-performance computing. MLegS is written in Modern Fortran; the code package is open-source under a BSD license.

Prior to its open-source release, MLegS was successfully used in several vortex dynamics studies in the context of wake vortices in the atmosphere. One example is S. Lee and P. S. Marcus (2023)[^2], where one can find the detailed mathematical formulation of the mapped Legendre (pseudo-)spectral method.

## Preparation
Refer to the Notebook tutorial to complete the prerequisites before using MLegS (`/tutorials/01_prerequisites.ipynb`). An online, ready-to-read (though non-interactive) version is available on the MLegS documentation site: [https://ucbcfd.github.io/MLegS/tutorial/prerequisites.html](https://ucbcfd.github.io/MLegS/tutorial/prerequisites.html).

## Links
- Code repository: [https://github.com/ucbcfd/MLegS](https://github.com/ucbcfd/MLegS)
- Documentation: [https://ucbcfd.github.io/MLegS](https://ucbcfd.github.io/MLegS)


## Contributors
- Jinge Wang (code - parallelization)
- Sangjoon Lee (code - modernization, documentation)
- [UC Berkeley CFD Lab](https://cfd.me.berkeley.edu) (codebase)

[^1]: [*J. Comput. Phys.*, vol. 137, no. 2, 1997, pp. 321-345](https://doi.org/10.1006/jcph.1997.5804)
[^2]: [*J. Fluid Mech.*, vol. 967, 2023, A2](https://doi.org/10.1017/jfm.2023.455)
