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
    """Extracts the prefix path for the given spack package"""

    if not use_bin:
        path = spec[package_name].prefix
    else:
        path = spec[package_name].prefix.bin
    path = path_replace(path, path_replacements)
    return path


def path_replace(path, path_replacements):
    """Replaces path key/value pairs from paht_replacements in path"""
    for key in path_replacements:
        path = path.replace(key,path_replacements[key])
    return path


class SeracLibraries(Package):
    """This is a set of libraries necessary for the developers of Serac.
       It also generates a host-config to be used for building Serac."""

    # TODO: Change this to a Spack BundledPackage that creates the host-config
    #       or change the git to be Serac's repo when it is on github

    #git = "ssh://git@cz-bitbucket.llnl.gov:7999/ser/serac.git"
    git = "https://github.com/LLNL/blt.git"

    version('develop', branch='develop', submodules=True, preferred=True)

    # List of Serac's third-party library dependencies
    depends_on('mfem +superlu-dist~shared')
    depends_on('superlu-dist~shared')
    depends_on('parmetis~shared')
    depends_on('metis~shared')

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
        compiler_string = str(spec.compiler).strip('%')
        host_config_filename = "{0}.cmake".format(compiler_string)
        host_config_path = os.path.abspath(os.path.join(env["SPACK_DEBUG_LOG_DIR"], host_config_filename))

        cfg = open(host_config_path, "w")
        cfg.write("####################################################################\n")
        cfg.write("# Generated host-config - Edit at own risk!\n")
        cfg.write("####################################################################\n")
        cfg.write("# Copyright (c) 2019, Lawrence Livermore National Security, LLC and\n")
        cfg.write("# other Serac Project Developers. See the top-level LICENSE file for\n")
        cfg.write("# details.\n")
        cfg.write("#\n")
        cfg.write("# SPDX-License-Identifier: (BSD-3-Clause) \n")
        cfg.write("####################################################################\n\n")

        cfg.write("#---------------------------------------\n")
        cfg.write("# SYS_TYPE: {0}\n".format(sys_type))
        cfg.write("# Compiler Spec: {0}\n".format(spec.compiler))
        cfg.write("# CMake executable path: %s\n" % cmake_exe)
        cfg.write("#---------------------------------------\n\n")

        #######################
        # Compiler Settings
        #######################

        cfg.write("#---------------------------------------\n")
        cfg.write("# Compilers\n")
        cfg.write("#---------------------------------------\n")
        cfg.write(cmake_cache_entry("CMAKE_C_COMPILER", c_compiler))
        cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER", cpp_compiler))

        #######################
        # MPI
        #######################

        cfg.write("#---------------------------------------\n")
        cfg.write("# MPI\n")
        cfg.write("#---------------------------------------\n")
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

        cfg.write("#---------------------------------------\n")
        cfg.write("# Library Dependencies\n")
        cfg.write("#---------------------------------------\n")

        path_replacements = {}

        # Try to find the common prefix of the TPL directory, including the compiler
        # If found, we will use this in the TPL paths
        compiler_str = str(spec.compiler).replace('@','-')
        prefix_paths = prefix.split( compiler_str )
        tpl_root = ""
        if len(prefix_paths) == 2:
            tpl_root = os.path.join( prefix_paths[0], compiler_str )
            path_replacements[tpl_root] = "${TPL_ROOT}"
            cfg.write(cmake_cache_entry("TPL_ROOT",tpl_root))

        mfem_dir = get_spec_path(spec, "mfem", path_replacements)
        cfg.write(cmake_cache_entry("MFEM_DIR", mfem_dir))

        #######################
        # Adding developer tools
        #######################

        #TODO: Change this to the common location (/usr/WS1/smithdev/tools) when thats available
        devtools_root = os.path.dirname(os.path.dirname(tpl_root))
        devtools_root = os.path.join(devtools_root, "devtools")

        if os.path.exists(devtools_root):
            cfg.write("#---------------------------------------\n")
            cfg.write("# Developer Tools\n")
            cfg.write("#---------------------------------------\n")

            #TODO: These paths should be able to come from a Spack upstream when they aren't tied to a spec
            astyle_path = "${DEVTOOLS_ROOT}/gcc-7.3.0/astyle-3.1/bin/astyle"
            cfg.write(cmake_cache_entry("ASTYLE_EXECUTABLE", astyle_path))

            cppcheck_path = "${DEVTOOLS_ROOT}/gcc-7.3.0/cppcheck-1.87/bin/cppcheck"
            cfg.write(cmake_cache_entry("CPPCHECK_EXECUTABLE", cppcheck_path))

            doxygen_path = "${DEVTOOLS_ROOT}/gcc-7.3.0/doxygen-1.8.15/bin/doxygen"
            cfg.write(cmake_cache_entry("DOXYGEN_EXECUTABLE", doxygen_path))

            sphinx_path = "${DEVTOOLS_ROOT}/gcc-7.3.0/py-sphinx-2.2.0/bin/sphinx-build"
            cfg.write(cmake_cache_entry("SPHINX_EXECUTABLE", sphinx_path))
        endif()

        #######################
        # Close and save
        #######################
        cfg.write("\n")
        cfg.close()

        # Fake install something so Spack doesn't complain
        mkdirp(prefix)
        install(host_config_path, prefix)
        print("Spack generated Serac host-config file: {0}".format(host_config_path))

