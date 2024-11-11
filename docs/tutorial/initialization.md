---
title: '#1. Initialization'
parent: Tutorials
nav_order: 2
---

# MLegS Tutorial 01: Initialization
*Disclaimer: This MLegS tutorial assumes a Linux or other Unix-based environment that supports bash terminal commands. If you are using Windows, consider installing the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install).*

In this tutorial, you will learn the initialization process for MLegS. Within a cylindrical coordinate system \\((r, \phi, z)\\), MLegS discretizes the radial direction \\((0 \le r < \infty)\\) using *mapped Legendre functions* ([Lee & Marcus, *J. Fluid Mech.*, 2023](https://doi.org/10.1017/jfm.2023.455)) as spectral basis elements. To achieve this, MLegS pre-calculates these functions with appropriate normalization at each collocation point, which is then pre-loaded into the global spectral transformation kit, `tfm`. This tutorial will guide you through the following steps:

1. **Compile and Run a 'Bare-Bones' MLegS Program Template**  
   - Compile and execute a program template located at `[root_dir]/src/apps/barebone_template.f90`, which only includes essential initialization steps.
2. **Experience MPI-Based Parallelism**
   - Although initializing the spectral transformation kit is a one-time process per program run, it can require substantial computational time as numerical resolution increases.
   - MLegS supports parallel computation for evaluating multiple spectral basis functions, enabling users to reduce computation time by leveraging MPI-based parallelism.
3. **Understand the Functional Approach of MLegS - What's Happening Behind the Scenes?**  
   - Gain a brief overview of the *mapped Legendre functions*, including how they serve as the spectral basis elements for constructing arbitrary scalar fields in three-dimensional space.

Completing this tutorial will help you understand the initialization steps in MLegS and the importance of this process for setting up simulations.

---

## Compile and Run a 'Bare-Bones' MLegS Program Template

With all the preliminary tasks from the previous tutorial successfully completed, you are now ready to compile a Fortran program that utilizes the functions and subroutines provided by the MLegS modules. According to the provided `Makefile` instructions, all program `.f90` files should be placed in the `[root_dir]/src/apps/` directory. The 'bare-bones' program template, `barebone_template.f90`, included in the package, is intended for advanced users who want to create custom MLegS programs. However, it is also helpful for any user to understand the essential initialization steps in MLegS before proceeding with iterative computations.

Run the following command to compile `barebone_template.f90` and generate an executable file, which will be saved in the `[root_dir]/build/bin/` directory.


```bash
cd ../ # Navigate to the root directory, assuming this notebook is opened in the default directory ([root_dir]/tutorials/).
```


```bash
# Do program compilation. You will be able to see what actually this Makefile command runs in the output. 
make barebone_template
```


```bash
# Once successfully compiled, the 'barebone_template' executable file must be seen:
ls -all ./build/bin/ | grep barebone_template
```

The compilation command should look like this:

`mpifort [option flags (-std08, etc.)] -J[mod_dir] -I[ext_inc_dir] -L[ext_lib_dir] -o [exe_file] [dependent_obj_files] -l[dependent_ext_lib_files]`

If you are using the Intel Fortran compiler, the syntax may vary slightly, but overall it should be similar. For details on compiler flags, refer to gfortran's [option summary](https://gcc.gnu.org/onlinedocs/gfortran/Option-Summary.html) (or ifx's [compiler options](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2025-0/compiler-options-001.html)). Advanced users can customize compiler options to optimize performance by modifying in `Makefile`, but note that options beyond those provided in the package have not been tested.

One key option flag is `-std08` (or `-std=f2008`), which ensures that the program and all dependent Fortran module files conform to the modern Fortran 2008 standard instead of the legacy F77 standard, which is one of the major updates made during the transition of MLegS from its original, proprietary in-house codebase to an open-source package.

Now let's run the program using single-core execution:


```bash
mpirun.openmpi -n 1 ./build/bin/barebone_template
# # for ifx + IntelMPI
# mpirun -n 1 ./build/bin/barebone_template
```

If the output appears as shown below, it means that the program has run successfully.

```
 Program Started
 @ XXXX-XX-XX XX:XX:XX
 
 Initialization
 elapsed time:        9.877282seconds
 
 *******************************************
 This program is solely a barebone template.
 No further computations shall be done.
 *******************************************
 
 Finalization
 elapsed time:        0.045591seconds
 
 Program Finished
 elapsed time:        0.000007seconds
 @ XXXX-XX-XX YY:YY:YY
```

---

## Experience MPI-Based Parallelism

As you may have noticed in the previous program run, the initialization process requires a significant amount of time due to the evaluation of numerous mapped Legendre functions across varying degrees and orders. The default input parameters provided in `[root_dir]/input.params` are \\( NRCHOP = 200 \\) and \\( NPCHOP = 65 \\), resulting in a total of \\( \sum_{m=0}^{65-1} (200 - m) = 10920 \\) mapped Legendre functions to be precomputed and loaded (details on this will be covered in the final section of this tutorial).

The good news is that these computations are *embarrassingly* parallelizable: each function can be evaluated independently with no depedence on the other functions' evaluations, allowing you to fully utilize your system's multi-core capability. In the following command, you’ll execute the program using dual-core processing, where each core handles half of the workload that a single core managed previously.


```bash
mpirun.openmpi -n 2 ./build/bin/barebone_template
# # for ifx + IntelMPI
# mpiexec -n 2 ./build/bin/barebone_template
```

Compared to the previous single-core run, you should notice that the computation time for the initialization stage has been nearly halved, indicating a nearly 2x speedup. The speedup is not a perfect 2x due to inter-processor communication overhead; each processor must share and transfer its locally computed evaluation data with the others so that all MPI workers have access to the fully pre-loaded mapped Legendre function information.

If your system has more than two processors, you can increase the number of MPI workers by adjusting `-n 2` to `-n [total # of processors]`.


```bash
# get the total number of processors of your system
np=$(nproc)
echo "The system's total number of processors is $np"
# run the program with your system's full multi-core capacity. 
#'--oversubscribe' is to make sure that the command below can be executed while bypassing the no. of available slots that may be limited by the default setup of OpenMPI.
# However, this flag is essentially unnecessary.
mpirun.openmpi -n $np --oversubscribe ./build/bin/barebone_template
# # for ifx + IntelMPI
# mpiexec -n $np ./build/bin/barebone_template
```

---

## Understand the Functional Approach of MLegS - What's Happening Behind the Scenes?

The associated Legendre functions are solutions to the general Legendre equation:

$$ \frac{d}{dx} \left[ (1-x^2) \frac{d}{dx} P_n^m (x) \right] + \left[ n(n+1) - \frac{m^2}{1-x^2} \right] P_n^m (x) = 0, $$

where the indices \\(n \\) and  \\(m \\) are the degree and order of the associated Legendre function \\( P_n^m (x) \\). In MLegS, these functions are algebraically mapped from \\( x \in [-1, 1) \\) to \\( r \in [0, \infty) \\) using the map parameter \\( L > 0 \\) such that \\( x = (r^2 + L^2)/(r^2 - L^2) \\). The mapped Legendre functions thus vary with the choice of \\( L \\) and are often denoted \\( P_{L_n}^m (r) \\). Now, \\( P_{L_n}^m (r) \\) is the solution to the following second-order differential equation:

$$ \frac{d}{dr} \left[ r \frac{d}{dr} P_{L_n}^m (r) \right] - \frac{m^2}{r} P_{L_n}^m (r) + \frac{4n(n+1)L^2 r}{(L^2 + r^2)^2} P_{L_n}^m (r) = 0. $$

According to Sturm-Liouville theory, the mapped Legendre functions \\( P_{L_{\|m\|}}^m \\), \\( P_{L_{\|m\|+1}}^m \\), \\( P_{L_{\|m\|+2}}^m \\), \\( \cdots \\) form a basis for any arbitrary function \\( f(r) \\) that decays sufficiently as \\( r \\) approaches infinity. For a polar function \\( f(r, \phi) \\) expressed as a Fourier series in \\( \phi \\), i.e., \\( f(r, \phi) = \sum_{m = -\infty}^{\infty} f_m (r) \exp(\rm{i}m \phi) \\), the following expansion in \\( r \\) is particularly useful:

$$ f(r, \phi) = \sum_{m = -\infty}^{\infty} \sum_{n = |m|}^{\infty} f_n^m P_{L_n}^m (r) \exp(\rm{i}m \phi), $$

as its truncated series inherently satisfies **analyticity** at the origin and decays quickly as \\( r \rightarrow \infty \\), regardless of where the truncation occurs, thus eliminating issues with coordinate singularities.

MLegS leverages this mapped Legendre expansion to discretize an arbitrary **real-valued** scalar field in a radially unbounded domain \\( 0 \le z < Z \\), e.g., \\( s (r, \phi , z) \\), as follows:

$$ s(r, \phi, z) = \sum_{\kappa = -NZCHOP+2}^{NZCHOP-1} \; \sum_{m = -NPCHOP+2}^{NPCHOP-1} \; \sum_{n = |m|}^{NRCHOP-1} f_n^{m\kappa} P_{L_n}^m (r) \exp(\rm{i}m \phi + \rm{i} \kappa z), $$

where \\( \kappa \equiv 2\pi / Z \\) represents the axial wavelength. The scalar field in the continuous domain \\( (r, \phi, z) \\) is now expressed as a finite set of discrete function coefficients \\( \left\\{ f_n^{m\kappa} \right\\} \\). Three truncation parameters—\\( NRCHOP \\), \\( NPCHOP \\), and \\( NZCHOP \\)—determine the accuracy of the discrete approximation; higher values yield greater accuracy. For efficient computation, MLegS restricts these parameters to meet the following conditions:
- \\( NRCHOP \\) must be even.
- \\( NPCHOP - 1 \\) must be even, with prime factors limited to 2, 3, and 5.
- \\( NZCHOP - 1 \\) must be even, with prime factors limited to 2, 3, and 5.

The discretized version of \\( s(r, \phi, z) \\) requires pre-loaded information in the form of \\( P_{L_n}^m (r) \\) for \\( m = 0, \cdots, NPCHOP-1 \\) and \\( n = \|m\|, \cdots, NRCHOP-1 \\) for each \\( m \\). The following table visually describes the required \\( P_{L_n}^m (r) \\) for a given \\( NRCHOP \\) and \\( NPCHOP \\):

|| Order 0 | Order 1 | Order 2 | \\( \cdots \\) | Order \\( NPCHOP-1 \\) |
|:--:|:--:|:--:|:--:|:--:|:--:|
|**Degree \\( NRCHOP-1 \\)**| \\( P_{L_{NRCHOP-1}}^{0} \\) | \\( P_{L_{NRCHOP-1}}^{1} \\) | \\( P_{L_{NRCHOP-1}}^{2} \\) | \\( \cdots \\) | \\( P_{L_{NRCHOP-1}}^{NPCHOP-1} \\) |
|**Degree \\( NRCHOP-2 \\)**| \\( P_{L_{NRCHOP-2}}^{0} \\) | \\( P_{L_{NRCHOP-2}}^{1} \\) | \\( P_{L_{NRCHOP-2}}^{2} \\) | \\( \cdots \\) | \\( P_{L_{NRCHOP-2}}^{NPCHOP-1} \\) |
|\\( \vdots \\)| \\( \vdots \\) | \\( \vdots \\) | \\( \vdots \\) | \\( \ddots \\) | \\( \vdots \\) |
|**Degree \\( NPCHOP-1 \\)**| \\( P_{L_{NPCHOP-1}}^{0} \\) | \\( P_{L_{NPCHOP-1}}^{1} \\) | \\( P_{L_{NPCHOP-1}}^{2} \\) | \\( \cdots \\) | \\( P_{L_{NPCHOP-1}}^{NPCHOP-1} \\) |
|\\( \vdots \\)| \\( \vdots \\) | \\( \vdots \\) | \\( \vdots \\) | \\( \cdots \\) | - |
|**Degree 2**| \\( P_{L_{2}}^{0} \\) | \\( P_{L_{2}}^{1} \\) | \\( P_{L_{2}}^{2} \\) | \\( \cdots \\) | - |
|**Degree 1**| \\( P_{L_{1}}^{0} \\) | \\( P_{L_{1}}^{1} \\) | - | \\( \cdots \\) | - |
|**Degree 0**| \\( P_{L_{0}}^{0} \\) | - | - | \\( \cdots \\) | - |
|Total number of elements| \\( NRCHOP \\) | \\( NRCHOP-1 \\) | \\( NRCHOP-2 \\) | \\( \cdots \\) | \\( NRCHOP-NPCHOP+1 \\) |

Accordingly, during the initialization stage of MLegS, a total of \\( \sum_{m=0}^{NPCHOP-1} (NRCHOP - m) \\) mapped Legendre functions must be evaluated.

To represent \\( s(r, \phi, z) \\) in its original physical form in the cylindrical domain \\( 0 \le z < Z \\), rather than its spectral coefficients, the physical domain is also discretized using the parameters \\( NR \\), \\( NP \\), and \\( NZ \\). The radial, azimuthal, and axial coordinates are discretized as:

$$ 0 < r_1 < r_2 < \cdots < r_{NR} < \infty, $$

$$ 0 = \phi_1 < \phi_2 < \cdots < \phi_{NP} < 2\pi, $$

$$ 0 = z_1 < z_2 < \cdots < z_{NZ} < Z, $$

where \\( r_i \\) values are collocated using the Gauss-Legendre quadrature rule ([Lee & Marcus, *J. Fluid Mech.*, 2023](https://doi.org/10.1017/jfm.2023.455); see also [Press *et al.*, *Numerical Recipes* (3rd edn.), 2007](https://numerical.recipes/) for its implementation), while \\( \phi_j \\) and \\( z_k \\) are equispaced. To prevent aliasing errors, these physical discretization parameters must satisfy:

- \\( NR \ge NRCHOP \\)
- \\( NP \ge 2(NPCHOP-1) \\)
- \\( NZ \ge 2(NZCHOP-1) \\)

For optimal use of numerical resources (avoiding both oversampling and aliasing), it is recommended to use the equality of each side of these inequalities.

---