## External packages

MLegS includes the following copies of external open-source code packages for necessary computing routines:

1. Multiple Precision Computation -- FM (v1.4) (Released 2024 Jan)
 - Source: [https://dmsmith.lmu.build/](https://dmsmith.lmu.build/)
 - Description: Fortran-based multiple-precision real, complex, and integer arithmetic
 - License: Refer to LICENSE in the directory (`./fmlib1.4/`)

2. Fast Fourier Transform Package -- FFTE (v7.0) (Released 2020 Aug)
 - Source: [http://www.ffte.jp/](http://www.ffte.jp/)
 - Description: Discrete Fourier Transforms of 1-, 2- and 3- dimensional sequences of length (2^p)(3^q)(5^r)
 - License: Refer to readme.txt in the directory (`./ffte-7.0/`)

3. Linear Algebra PACKage -- LAPACK (v3.12.0) (Released 2023 Nov)
 - Source: [https://www.netlib.org/lapack/](https://www.netlib.org/lapack/)
 - Description: Solving systems of linear equations, eigenvalue problems, and singular value problems
 - License: Refer to LICENSE in the directory (`./lapack/`)

4. General-purpose 2D pencil decomposition module -- 2DECOMP&FFT (v2.0.3) (Released 2024 Mar)
 - Source: [https://2decomp-fft.github.io/](https://2decomp-fft.github.io/)
 - Description: General-purpose 2D pencil decomposition for data distribution
 - License: Refer to LICENSE in the directory (`./2decomp-fft/`)
 - *Note: We do not make use of 2DECOMP&FFT's built-in FFT module; FFTE is our FFT engine.*

Please refer to the license file in each directory for the full open-source license text of these packages.

### How to build

As general Fortran is known to be non-ABI compatible, it is recommended to use the same Fortran compiler (as well as the corresponding C compiler, if necessary) that you intend to use to build and install MLegS when building these external libraries. Currently, GNU's gfortran (v11 or later) with OpenMPI (v4.1) and Intel's ifx (v2024 or later) with its companion MPI library have been tested and verified.

In the provided bash shell script (`CMake_build.sh`), change the following two arguments in the script if you want to switch the compiler; the default option is GNU.

```bash
...
FC="gfortran" # change gfortran to ifx if one wants to use Intel's oneAPI compiler
CC="gcc" # change gcc to icx if one wants to use Intel's oneAPI compiler
...
```

Once the compiler arguments are declared, simply run `./CMake_build.sh`. The whole build process is automated by [CMake](https://cmake.org/); if CMake is not installed on your machine, you may need to install it first.

After all the external packages are successfully built, the static library files will be created in `./lib/`, and the binary module files will be stored in `./inc/`. Note that both of these subdirectories are automatically created during the build process, so they don't need to be created in advance.