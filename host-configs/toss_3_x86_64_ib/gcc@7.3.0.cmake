# Copyright (c) 2019, Lawrence Livermore National Security, LLC and
# other Serac Project Developers. See the top-level LICENSE file for
# details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

# SYS_TYPE: toss_3_x86_64_ib
# Compiler Spec: gcc@7.3.0
##################################

# CMake executable path: /usr/tce/packages/cmake/cmake-3.9.2/bin/cmake


##############
# Developer Tools
##############
set(SERAC_DEVTOOLS_ROOT "/usr/WS2/white238/serac/repo/devtools/gcc-7.3.0" CACHE PATH "")

set(ASTYLE_EXECUTABLE "${SERAC_DEVTOOLS_ROOT}/astyle-3.1/bin/astyle")
set(CPPCHECK_EXECUTABLE "${SERAC_DEVTOOLS_ROOT}/cppcheck-1.87/bin/cppcheck")
set(DOXYGEN_EXECUTABLE "${SERAC_DEVTOOLS_ROOT}/doxygen-1.8.15/bin/doxygen")
set(SPHINX_EXECUTABLE "${SERAC_DEVTOOLS_ROOT}/py-sphinx-2.2.0/bin/sphinx-build")


##############
# Compilers
##############

# C compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-7.3.0/bin/gcc" CACHE PATH "")

# C++ compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-7.3.0/bin/g++" CACHE PATH "")

##############
# MPI
##############

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-7.3.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-7.3.0/bin/mpicxx" CACHE PATH "")

set(MPIEXEC "/usr/bin/srun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############
# Other machine specifics
##############

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")


