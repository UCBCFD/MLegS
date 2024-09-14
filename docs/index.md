---
title: Home
layout: home
---

# **MLegS**: Modernized and Parallelized **M**apped **Leg**endre **S**pectral Method Code
MLegS (**M**apped **Leg**endre **S**pectral Method Code) is a code package based on a parallelized spectral method for vortical (swirling) flow simulations in an unbounded domain.

Based on the numerical algorithm proposed by T. Matsushima and P. S. Marcus (1997)[^1], MLegS incorporates scalable multiprocessing interfaces for high-performance computing. MLegS is written in Modern Fortran; the code package is open-source under an MIT license.

Prior to its open-source release, MLegS was successfully used in several vortex dynamics studies in the context of wake vortices in the atmosphere. One example is S. Lee and P. S. Marcus (2023)[^2], where one can find the detailed mathematical formulation of the mapped Legendre (pseudo-)spectral method.

## Links
- Code repository: [https://github.com/ucbcfd/MLegS](https://github.com/ucbcfd/MLegS)
- Documentation: [https://ucbcfd.github.io/MLegS](https://ucbcfd.github.io/MLegS)


## Contributors
- Jinge Wang
- Sangjoon Lee
- [UC Berkeley CFD Lab](https://cfd.me.berkeley.edu)

[^1]: [*J. Comput. Phys.*, vol. 137, no. 2, 1997, pp. 321-345](https://doi.org/10.1006/jcph.1997.5804)
[^2]: [*J. Fluid Mech.*, vol. 967, 2023, A2](https://doi.org/10.1017/jfm.2023.455)