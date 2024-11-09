---
title: '#0. Prerequisites'
parent: Tutorials
nav_order: 1
---

# MLegS Tutorial #0. Prerequisites

In this tutorial, you’ll learn how to prepare your environment to work with the MLegS package. By completing these prerequisites, you will be ready to compile and run main programs that use MLegS. The following steps will be covered:
1. Configure a compiler and MPI (Message Passing Interface)
  - Ensure you have a compatible compiler installed (e.g., GNU or Intel)
  - Set up MPI (Message Passing Interface) on your machine, essential for parallel execution.
2. Set up external libraries required to run MLegS
  - MLegS relies on certain external libraries for matrix operations, fast Fourier transforms and multi-precision arithmetic.
  - Install necessary libraries such as LAPACK, FFTE and FM. You do not need to download them; these external libraries, with the compatibility checked for MLegS, are provided in the `[root_dir]/external/` directory.
3. Compile MLegS modules.
  - With your compiler, MPI, and external libraries set up, proceed to compile the core MLegS modules.
  - Use the provided `Makefile` instructions to compile the modules via one-line commands (e.g., `make mods`)

Completing this setup ensures that you’re fully prepared to work with programs that use the MLegS package.

---

## Configure a Compiler and MPI Interfaces

To ensure compatibility with MLegS, use one of the compiler and MPI setups. The following combinations are recommended for compiling and running MLegS:
1. GNU Fortran compiler `gfortran` <v11.2 or later> with OpenMPI <v4.1 or later>
2. Intel OneAPI Fortran compiler `ifx` <v2024.1.0 or later> with accompanying Intel MPI

While Intel's legacy Fortran compiler, `ifort`, may be used as an alternative, please note that full compatibility has not been tested or guaranteed.

You can install `gfortran` and OpenMPI via your system's package manager (`apt` for Ubuntu/Debian, `yum` for CentOS), as these packages are generally available in official repositories. 

For Intel compilers, download and install the toolkit from [Intel's official download page](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.hbhvru), selecting the distribution that matches your system. After installation, you need to source the Intel environment script (e.g., `source /opt/intel/oneapi/setvars.sh` to set environment variables.

Aftr installation, verify the versions with:


```bash
# for gfortran + OpenMPI
gfortran --version
mpirun.openmpi --version
# # for ifx + IntelMPI
# ifx --version
# mpirun --version
```

---

## Set up external libraries required to run MLegS

MLegS depends on several external open-source libraries that provide essential routines for computations. The following packages are required:
1. Multiple Precision Computation -- FM (v1.4) ([Source](https://dmsmith.lmu.build/))
2. Fast Fourier Transform Package -- FFTE (v7.0) ([Source](http://www.ffte.jp/))
3. Linear Algebra PACKage -- LAPACK (v3.12.0) ([Source](https://www.netlib.org/lapack/))

All required libraries are included in the `[root_dir]/external/` directory of the MLegS package, so manual downloading is not needed. MLegS provides a bash script to compile these libraries in one step, simplifying the setup process. Before running the compilation script, ensure the following prerequisites are met:
  - `CMake` <v3.10 or later> is installed on your system. You can check if `CMake` is installed and verify the version with:


```bash
cmake --version
```

- If it is not installed or outdated, use your package manager or visit `CMake`'s [official website](https://cmake.org/) to install or update it.

Once CMake is installed, navigate to the `[root_dir]/external/` directory and check if the provided compliation script (`CMake_build.sh`) is executable.


```bash
cd ../external/ # Navigate to the external directory, assuming this notebook is opened in the default directory ([root_dir]/tutorials/).
```


```bash
chmod +x ./CMake_build.sh # Add the executable option
# Now let's check that CMake_build.sh is actually executable. If you get something like '-rwxrwxr-x' of 3 x's, then it's executable.
ls -all | grep CMake_build.sh 
```

By default, the script uses `gfortran` (with its compatible C compiler) for compilation, but you can easily change the compiler by modifying the compiler specification within the script.


```bash
# in [root_dir]/external/CMake_build.sh, using any preferred text editor, update the following lines if one wants to use Intel
# ...
# FC = "gfortran" (-> replace it with "ifx", Intel's oneAPI Fortran compiler)
# CC = "gcc" (-> replace it with "icx", Intel's oneAPI C compiler
# ...
# or simply run this bash command:
#   sed -i 's/"gfortran"/"ifx"/g' ./CMake_build.sh
#   sed -i 's/"gcc"/"icx"/g' ./CMake_build.sh
# Now let's check what compilers are in use:
head -n 20 ./CMake_build.sh | tail -4
```

You are now ready to run the automated external library compilation script. Execute the following command:


```bash
# Since this script run yields very lengthy outputs, we suppress them here. 
# Generally this automated compilation task takes up to 5-10 minutes. Take a cup of coffee... :)
./CMake_build.sh > /dev/null 2>&1
```

If all external libraries are compiled successfully, you will see two new directories created in `[root_dir]/external/`: `./inc/` and `./lib/`. In `./lib/`, you should find the following **four** static library files: `libblas.a` and `liblapack.a` from LAPACK, `libfm.a` from FM, and `libffte.a` from FFTE. In `./inc/`, you should see several FM-related module files with a `.mod` extension (e.g., `fmvals.mod`, `fmzm.mod`, etc.).


```bash
# See if the library files are all generated and stored in the lib directory.
ls ./lib/ ./inc/
```

---

## Compile MLegS Modules

Now that the Fortran compiler, MPI, and the external libraries are set up, it is possible to proceed to compile the core MLegS modules. By default, the whole compilation process of MLegS is assumed to be done in its root directory. In the root directory, `Makefile` provides the instructions to compile the modules via a single line command, `make mods`.

Before doing so, move to the root directory and ensure that the `Makefile` instructions runs over the same Fortran compiler as what you have utilized for the external library compilation. By default, the instructions use `gfortran`.


```bash
cd ../ # Navigate to the root directory, assuming that you were at [root_dir]/external/.
```


```bash
# in [root_dir]/Makefile, using any preferred text editor, update the following lines if one wants to use Intel
# ...
# OMPI_FC = "gfortran" (-> replace it with "ifx", Intel's oneAPI Fortran compiler)
# ...
# or simply run this bash command:
#   sed -i 's/OMPI_FC = gfortran/OMPI_FC = ifx/g' ./Makefile
# Now let's check what compiler is in use:
head -n 4 ./Makefile | tail -4
```

MLegS, in its most recent version, includes a stack of seven modules designed for MPI-parallelized simulations within a radially unbounded computational domain. These module files are stored in `[root_dir]/src/modules/`, where you can find the header information for all functions and subroutines. The actual numerical calculations and I/O operations are implemented in submodule files located in `[root_dir]/src/submodules/`. While we won’t go into detail about each module's functionality in this tutorial, here’s a brief overview of what each module contains:
- mlegs_envir.f90: Defines environmental variables (e.g., default precision for integers, complex numbers, and real numbers, as well as simulation time information) and MPI parameters, along with basic auxiliary functions.
- mlegs_base.f90: Sets up global simulation parameters, including spectral element counts in each direction, time-stepping details, viscosity, and field data I/O.
- mlegs_misc.f90: Provides miscellaneous utilities, such as timers and generic Fortran array I/O functions.
- mlegs_genmat.f90: Contains routines for general matrix operations.
- mlegs_bndmat.f90: Includes operations for banded matrices (e.g., diagonal, tridiagonal).
- mlegs_spectfm.f90: Contains tools for spectral transformations and related spectral operations.
- mlegs_scalar.f90: Defines a distributed scalar class with spatially discretized operations and temporal advancement schemes.

Inter-module dependencies are specified in `Makefile.dep`, which the `Makefile` instructions use to determine the correct compilation order automatically This setup allows you to compile without worrying about module dependencies. Simply execute the following command: 


```bash
# 'clean' subcommand precedes in case there are old compiled files; this is not required and running only 'make mods' is generally sufficient.
make clean && make mods
```

If all modules compiled successfully, a new directory `[root_dir]/build/` will be generated. Inside it, you should find two subdirectories: `./mod/` for module files with `.mod` and `.smod` extensions, and `./obj/` with a `.o` extension.


```bash
# See if the library files are all generated and stored in the build directory.
ls ./build/mod/ ./build/obj/
```

With this, all preliminary tasks required to use MLegS are complete!

---