---
title: Home
layout: home
nav_order: 1
---

# **MLegS**: Modernized and Parallelized **M**apped **Leg**endre **S**pectral Method Code
MLegS (**M**apped **Leg**endre **S**pectral Method Code) is a code package based on a modernized and parallelized spectral method for numerical simulations in a radially unbounded domain.

Based on the numerical algorithm proposed by T. Matsushima and P. S. Marcus (1997)[^1], MLegS incorporates scalable multiprocessing interfaces for high-performance computing. MLegS is written in Modern Fortran; the code package is open-source under a BSD license.

Prior to its open-source release, MLegS was successfully used in several vortex dynamics studies in the context of wake vortices in the atmosphere. One example is S. Lee and P. S. Marcus (2023)[^2], where one can find the detailed mathematical formulation of the mapped Legendre (pseudo-)spectral method.

---

## Contributors
- Jinge Wang (code - parallelization)
- Sangjoon Lee (code - modernization, documentation)
- [UC Berkeley CFD Lab](https://cfd.me.berkeley.edu) (codebase)

---

[^1]: Matsushima, T., & Marcus, P. S. (1997). A spectral method for unbounded domains. Journal of Computational Physics, 137(2), 321–345. [https://doi.org/10.1006/jcph.1997.5804](https://doi.org/10.1006/jcph.1997.5804)
[^2]: Lee, S., & Marcus, P. S. (2023). Linear stability analysis of wake vortices by a spectral method using mapped Legendre functions. Journal of Fluid Mechanics, 967, A2. [https://doi.org/10.1017/jfm.2023.455](https://doi.org/10.1017/jfm.2023.455)