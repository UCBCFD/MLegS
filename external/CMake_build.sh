#!/bin/bash
# To run this bash script, make sure that cmake is installed in your system.
# Using apt, you may install cmake by running 'apt install cmake'
# (re)create a lib directory where all external libraries' symbolic links are
if [ -d ./lib ]; then
	rm -rf ./lib
fi
mkdir -p ./lib
if [ -d ./inc ]; then
	rm -rf ./inc
fi
mkdir -p ./inc

# Some packages require both the Fortran and C compilers, requiring interoperability to each other.
# Thus, it is recommended to explicitly specify what compilers you want to use in case your system
# possesses multiple compilers from different distributors (e.g., GNU and Intel).
# Change the COMPILER pair e.g., FC=ifx && CC=icx, FC=gfortran && CC=gcc, etc.
FC="gfortran" # change gfortran to ifx if one wants to use Intel's oneAPI compiler, and vice versa
CC="gcc" # change gcc to icx if one wants to use Intel's oneAPI compiler, and vice versa


###  generate the ffte library
cd ./ffte-7.0
# (re)create and go into a new build directory
if [ -d ./build ]; then
	rm -rf ./build
fi
mkdir -p ./build
cd ./build

# run cmake to create a script to build a library
cmake -DCMAKE_Fortran_COMPILER=${FC} ..

# run the build
cmake --build . -j

ln -s $(pwd)/lib/libffte.a ../../lib/libffte.a
cd ../../

### generate the fm library
cd ./fm1.4
# (re)create and go into a new build directory
if [ -d ./build ]; then
	rm -rf ./build
fi
mkdir -p ./build
cd ./build

# run cmake to create a script to build a library
cmake -DCMAKE_Fortran_COMPILER=${FC} ..

# run the build
cmake --build . -j

ln -s $(pwd)/lib/libfm.a ../../lib/libfm.a
ln -s $(pwd)/inc/* ../../inc/
cd ../../

### generate the lapack library
# Note: As lapack is included as a submodule, it's not downloaded when using the github release.
if [ ! -d "./lapack" ] || [ -z "$(ls -A "./lapack" 2>/dev/null)" ]; then
    rm -rf ./lapack
    git clone --depth 1 --branch v3.12.0 https://github.com/Reference-LAPACK/lapack.git
fi
# The above commands help download the lapack module (v3.12.0) from its official repo, in case
# ./lapack does not exist or even if it exists, the directory is empty (i.e., not downloaded).
cd ./lapack
# (re)create and go into a new build directory
if [ -d ./build ]; then
	rm -rf ./build
fi
mkdir -p ./build
cd ./build

# run cmake to create a script to build a library
cmake -DCMAKE_Fortran_COMPILER=${FC} -DCMAKE_C_COMPILER=${CC} -DCMAKE_INSTALL_LIBDIR=$(pwd)/lib \
      -DBUILD_SINGLE=OFF -DBUILD_DOUBLE=ON -DBUILD_COMPLEX=OFF -DBUILD_COMPLEX16=ON ..

# run the build
cmake --build . -j --target install

ln -s $(pwd)/lib/libblas.a ../../lib/libblas.a
ln -s $(pwd)/lib/liblapack.a ../../lib/liblapack.a
cd ../../

exit 0
