# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install serac
#
# You can edit this file again by typing:
#
#     spack edit serac
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *

import socket
import os

import llnl.util.tty as tty
from os import environ as env


def cmake_cache_entry(name, value, comment=""):
    """Generate a string for a cmake cache variable"""

    return 'set(%s "%s" CACHE PATH "%s")\n\n' % (name,value,comment)


def cmake_cache_option(name, boolean_value, comment=""):
    """Generate a string for a cmake configuration option"""

    value = "ON" if boolean_value else "OFF"
    return 'set(%s %s CACHE BOOL "%s")\n\n' % (name,value,comment)


def get_spec_path(spec, package_name, path_replacements = {}, use_bin = False) :
    """Extracts the prefix path for the given spack package
       path_replacements is a dictionary with string replacements for the path.
    """

    if not use_bin:
        path = spec[package_name].prefix
    else:
        path = spec[package_name].prefix.bin

    for key in path_replacements:
        path = path.replace(key,path_replacements[key])

    return path

class SeracLibraries(Package):
    """This is a set of libraries necessary for the developers of Serac.
       It also generates a host-config to be used for building Serac."""
    #git = "ssh://git@cz-bitbucket.llnl.gov:7999/ser/serac.git"
    git = "https://github.com/LLNL/blt.git"

    version('develop', branch='develop', submodules=True, preferred=True)

    # List of Serac's third-party library dependencies
    depends_on('mfem +superlu-dist')


    def install(self, spec, prefix):

        #######################
        # Compiler Info
        #######################
        c_compiler = env["SPACK_CC"]
        cpp_compiler = env["SPACK_CXX"]


        #######################################################################
        # By directly fetching the names of the actual compilers we appear
        # to doing something evil here, but this is necessary to create a
        # 'host config' file that works outside of the spack install env.
        #######################################################################

        sys_type = spec.architecture
        # if on llnl systems, we can use the SYS_TYPE
        if "SYS_TYPE" in env:
            sys_type = env["SYS_TYPE"]

        ##############################################
        # Find and record what CMake is used
        ##############################################

        cmake_exe = spec['cmake'].command.path
        compiler_string = str(spec.compiler).strip('%').replace('@','_').replace('.','_')
        host_cfg_fname = "%s__%s__serac.cmake" % (sys_type,
                                                  compiler_string)

        cfg = open(host_cfg_fname, "w")
        cfg.write("##################################\n")
        cfg.write("# Spack generated host-config\n")
        cfg.write("##################################\n")
        cfg.write("# {0}-{1}\n".format(sys_type, spec.compiler))
        cfg.write("##################################\n\n")

        # Include path to cmake for reference
        cfg.write("# cmake executable path: %s\n\n" % cmake_exe)

        #######################
        # Compiler Settings
        #######################

        cfg.write("#######\n")
        cfg.write("# Using %s compiler spec\n" % spec.compiler)
        cfg.write("#######\n\n")
        cfg.write(cmake_cache_entry("CMAKE_C_COMPILER", c_compiler))
        cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER", cpp_compiler))

        #######################
        # MPI
        #######################

        cfg.write(cmake_cache_entry("ENABLE_MPI", "ON"))
        cfg.write(cmake_cache_entry("MPI_C_COMPILER", spec['mpi'].mpicc))
        cfg.write(cmake_cache_entry("MPI_CXX_COMPILER",
                                    spec['mpi'].mpicxx))
        mpiexe_bin = join_path(spec['mpi'].prefix.bin, 'mpiexec')
        if os.path.isfile(mpiexe_bin):
            # starting with cmake 3.10, FindMPI expects MPIEXEC_EXECUTABLE
            # vs the older versions which expect MPIEXEC
            if self.spec["cmake"].satisfies('@3.10:'):
                cfg.write(cmake_cache_entry("MPIEXEC_EXECUTABLE",
                                            mpiexe_bin))
            else:
                cfg.write(cmake_cache_entry("MPIEXEC",
                                                mpiexe_bin))

        #######################
        # Adding dependencies
        #######################

        path_replacements = {}

        # Try to find the common prefix of the TPL directory, including the compiler
        # If found, we will use this in the TPL paths
        compiler_str = str(spec.compiler).replace('@','-')
        prefix_paths = prefix.split( compiler_str )
        if len(prefix_paths) == 2:
            tpl_root = pjoin( prefix_paths[0], compiler_str )
            path_replacements[tpl_root] = "${TPL_ROOT}"
            cfg.write("# Root directory for generated TPLs\n")
            cfg.write(cmake_cache_entry("TPL_ROOT",tpl_root))

        mfem_dir = get_spec_path(spec, "mfem", path_replacements)
        cfg.write(cmake_cache_entry("MFEM_DIR", mfem_dir))

        #######################
        # Adding developer tools
        #######################

        cfg.write("#######\n")
        cfg.write("# Developer tools\n")
        cfg.write("#######\n\n")

        #TODO: this should be in a common location
        devtools_root = "/usr/WS2/white238/serac/repo/devtools/gcc-7.3.0"
        path_replacements[devtools_root] = "${DEVTOOLS_ROOT}"
        cfg.write(cmake_cache_entry("DEVTOOLS_ROOT", devtools_root))

        #TODO: Can this come from a upstream
        cfg.write(cmake_cache_entry("ASTYLE_EXECUTABLE", "/usr/WS2/white238/serac/repo/devtools/gcc-7.3.0/astyle-3.1/bin/astyle"))
        cfg.write(cmake_cache_entry("CPPCHECK_EXECUTABLE", "/usr/WS2/white238/serac/repo/devtools/gcc-7.3.0/cppcheck-1.87/bin/cppcheck"))
        cfg.write(cmake_cache_entry("DOXYGEN_EXECUTABLE", "/usr/WS2/white238/serac/repo/devtools/gcc-7.3.0/doxygen-1.8.15/bin/doxygen"))
        cfg.write(cmake_cache_entry("SPHINX_EXECUTABLE", "/usr/WS2/white238/serac/repo/devtools/gcc-7.3.0/py-sphinx-2.2.0/bin/sphinx-build"))


        #######################
        # Close host-config
        #######################

        cfg.write("##################################\n")
        cfg.write("# End Spack generated host-config\n")
        cfg.write("##################################\n")
        cfg.close()

        host_cfg_fname = os.path.abspath(host_cfg_fname)
        tty.info("Spack generated Serac host-config file: " + host_cfg_fname)

