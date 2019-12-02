# Copyright (c) 2019, Lawrence Livermore National Security, LLC and
# other Serac Project Developers. See the top-level LICENSE file for
# details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

# SYS_TYPE: toss_3_x86_64_ib
# Compiler Spec: clang@4.0.0
##################################

# CMake executable path: /usr/tce/packages/cmake/cmake-3.8.2/bin/cmake

##############
# Compilers
##############

# C compiler used by spack
set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-4.0.0/bin/clang" CACHE PATH "")

# C++ compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-4.0.0/bin/clang++" CACHE PATH "")

##############
# MPI
##############

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.2-clang-4.0.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.2-clang-4.0.0/bin/mpicxx" CACHE PATH "")

set(MPIEXEC "/usr/bin/srun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############
# Other machine specifics
##############

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")


