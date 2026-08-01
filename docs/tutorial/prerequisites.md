---
title: '#1. Prerequisites'
parent: Tutorials
nav_order: 1
---

# MLegS Tutorial 01: Prerequisites
*Disclaimer: This MLegS tutorial assumes a Unix-based environment with a bash-compatible terminal. Linux and macOS are supported; on Windows, consider installing the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install).*

In this tutorial, you’ll learn how to prepare your environment to work with the MLegS package. By completing these prerequisites, you will be ready to compile and run main programs that use MLegS. The following steps will be covered:

1. **Configure a compiler and MPI (Message Passing Interface)**
   - Ensure you have a compatible compiler installed (e.g., GNU or Intel)
   - Set up MPI (Message Passing Interface) on your machine, essential for parallel execution.
2. **Set up external libraries required to run MLegS**
   - MLegS relies on certain external libraries for matrix operations, fast Fourier transforms and multi-precision arithmetic.
   - Install necessary libraries such as LAPACK, FFTE and FM. You do not need to download them; these external libraries, with the compatibility checked for MLegS, are provided in the `[root_dir]/external/` directory. 
3. **Compile MLegS modules**
   - With your compiler, MPI, and external libraries set up, proceed to compile the core MLegS modules.
   - Use the provided `Makefile` instructions to compile the modules via one-line commands (e.g., `make mods`)


Completing this tutorial ensures that you’re fully prepared to work with programs that use the MLegS package.

---

## Table of Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---


## Configure a Compiler and MPI Interfaces

To ensure compatibility with MLegS, use one of the compiler and MPI setups. The following combinations are recommended for compiling and running MLegS:
1. GNU Fortran compiler `gfortran` <v11.2 or later> with OpenMPI <v4.1 or later>
2. Intel OneAPI Fortran compiler `ifx` <v2024.1.0 or later> with accompanying Intel MPI

While Intel's legacy Fortran compiler, `ifort`, may be used as an alternative, please note that full compatibility has not been tested or guaranteed.

You can install `gfortran` and OpenMPI via your system's package manager (`apt` for Ubuntu/Debian, `yum` for CentOS, or Homebrew on macOS):

```bash
# macOS (Homebrew)
brew install gcc open-mpi cmake
```

For Intel compilers, download and install the toolkit from [Intel's official download page](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.hbhvru), selecting the distribution that matches your system. After installation, you need to source the Intel environment script (e.g., `source /opt/intel/oneapi/setvars.sh`) to set environment variables.

After installation, verify the versions with:


```bash
#! bash
# for gfortran + OpenMPI
gfortran --version
mpirun --version
# # for ifx + IntelMPI
# source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1; export PATH="/opt/intel/oneapi:$PATH"
# ifx --version
# mpiexec --version
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
#! bash
cmake --version
```

- If it is not installed or outdated, use your package manager or visit `CMake`'s [official website](https://cmake.org/) to install or update it.

Once CMake is installed, navigate to the `[root_dir]/external/` directory and check if the provided compliation script (`CMake_build.sh`) is executable.


```bash
#! bash
cd external/ # Navigate to the external directory from the repository root.
```


```bash
#! bash
chmod +x ./CMake_build.sh # Add the executable option
# Now let's check that CMake_build.sh is actually executable. If you get something like '-rwxrwxr-x' of 3 x's, then it's executable.
ls -all | grep CMake_build.sh 
```

By default, the script uses `gfortran` and a compatible C compiler.  The
compiler and parallel build count can be overridden without editing the
script; this is useful on macOS, where Homebrew provides GNU Fortran and
Apple Clang as `cc`.


```bash
#! bash
# macOS with Homebrew GNU Fortran:
FC=gfortran CC=cc BUILD_JOBS=2 ./CMake_build.sh

# Intel oneAPI (on a supported Linux installation):
# FC=ifx CC=icx BUILD_JOBS=2 ./CMake_build.sh
```

You are now ready to run the automated external library compilation script. Execute the following command:


```bash
#! bash
# The script now stops at the first failed dependency build.  Keep the output
# visible while checking a new compiler installation.
./CMake_build.sh
```

If all external libraries are compiled successfully, you will see two new directories created in `[root_dir]/external/`: `./inc/` and `./lib/`. In `./lib/`, you should find the following **four** static library files: `libblas.a` and `liblapack.a` from LAPACK, `libfm.a` from FM, and `libffte.a` from FFTE. In `./inc/`, you should see several FM-related module files with a `.mod` extension (e.g., `fmvals.mod`, `fmzm.mod`, etc.).


```bash
#! bash
# See if the library files are all generated and stored in the lib directory.
ls ./lib/ ./inc/
```

---

## Compile MLegS Modules

Now that the Fortran compiler, MPI, and the external libraries are set up, it is possible to proceed to compile the core MLegS modules. By default, the whole compilation process of MLegS is assumed to be done in its root directory. In the root directory, `Makefile` provides the instructions to compile the modules via a single line command, `make mods`.

Before doing so, move to the root directory and ensure that the `Makefile`
uses the same Fortran compiler as the external-library build.  GNU Fortran is
the default; the compiler wrapper can be selected with `OMPI_FC`.


```bash
#! bash
if [ -f Makefile ]; then ROOT_DIR="$PWD"; else ROOT_DIR="$(cd .. && pwd)"; fi
cd "$ROOT_DIR"
```


```bash
#! bash
# GNU Fortran (default, including macOS/Homebrew):
make mods

# Intel oneAPI on a supported Linux installation:
# make OMPI_FC=ifx mods
```

The GNU Fortran Makefile path does not pass Linux-only `-mcmodel=medium`
flags, so the same command works with Homebrew GCC on Apple silicon and Intel
macOS. Site-specific flags can be appended with
`GFORTRAN_EXTRA_FFLAGS="..."`.

MLegS, in its most recent version, includes a stack of seven modules designed for MPI-parallelized simulations within a radially unbounded computational domain. These module files are stored in `[root_dir]/src/modules/`, where you can find the header information for all functions and subroutines. The actual numerical calculations and I/O operations are implemented in submodule files located in `[root_dir]/src/submodules/`. While we won’t go into detail about each module's functionality in this tutorial, here’s a brief overview of what each module contains:

| Module name | Description|
|:-- |:-- |
| mlegs_envir.f90 | Defines environmental variables (e.g., default precision for integers, complex numbers, and real numbers, as well as simulation time information) and MPI parameters, along with basic auxiliary functions. |
| mlegs_base.f90 | Sets up global simulation parameters, including spectral element counts in each direction, time-stepping details, viscosity, and field data I/O. |
| mlegs_misc.f90 | Provides miscellaneous utilities, such as timers and generic Fortran array I/O functions. |
| mlegs_genmat.f90 | Contains routines for general matrix operations. |
| mlegs_bndmat.f90 | Includes operations for banded matrices (e.g., diagonal, tridiagonal). |
| mlegs_spectfm.f90 | Contains tools for spectral transformations and related spectral operations. |
| mlegs_scalar.f90 | Defines a distributed scalar class with spatially discretized operations and temporal advancement schemes. |

Inter-module dependencies are specified in `Makefile.dep`, which the `Makefile` instructions use to determine the correct compilation order automatically This setup allows you to compile without worrying about module dependencies. Simply execute the following command: 


```bash
#! bash
# 'clean' subcommand precedes in case there are old compiled files; this is not required and running only 'make mods' is generally sufficient.
make clean && make mods
```

If all modules compiled successfully, a new directory `[root_dir]/build/` will be generated. Inside it, you should find two subdirectories: `./mod/` for module files with `.mod` and `.smod` extensions, and `./obj/` with a `.o` extension.


```bash
#! bash
# See if the library files are all generated and stored in the build directory.
ls ./build/mod/ ./build/obj/
```

With this, all preliminary tasks required to use MLegS are complete!

## Validate all tutorials

Before submitting a change, run the repository tutorial gate from the root directory. It checks the documentation and notebook command surfaces, builds every documented example, and runs the tutorial chains in an isolated temporary directory:

```bash
python3 tools/validate_tutorials.py
```

The same gate runs automatically for pushes and pull requests in GitHub Actions. Use `--skip-build` when the current `build/bin/` executables are already available. Full notebook-cell execution is available with `--execute-notebooks` when `nbformat`, `nbclient`, and the plotting dependencies are installed.
