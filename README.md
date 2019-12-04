Serac
====

Serac is a 3D implicit nonlinear thermal-structural simulation code. It's primary purpose is to investigate multiphysics abstraction strategies and implicit finite element-based alogrithm development for emerging computing architectures. It also serves as a proxy-app for LLNL's DIABLO and ALE3D codes.

Getting Started
------

Serac uses git submodules, to clone the project:

```
git clone --recursive ssh://git@cz-bitbucket.llnl.gov:7999/ser/serac.git
```

The easiest path to install both Serac and its dependencies is to use Spack. This has been encapsulated
using Uberenv. This will download and configure Spack and build all of Serac's dependencies. It also
generates a host-config file that has all the necessary build information for Serac. The CMake
configuration phase has also been encapsulated in config-build.py.

If you want the optional developer tools, run the following command:

```
`scripts/uberenv/uberenv.py --package=serac_devtools --prefix=`pwd`/../devtools --spec=%gcc@8.1.0`
```

This will enable Doxygen, Sphinx, Cppcheck, and AStyle in your build.  This is not required.

Here is an example of how to build Serac and its build dependencies with clang 4.0.0:

1. `scripts/uberenv/uberenv.py --prefix=`pwd`/../libs --spec=%clang@4.0.0`
2. `./config-build.py -hc <library build directory>/clang@4.0.0.cmake`
3. `cd build-<system-and-toolchain>`
4. `cmake --build .`
5. `ctest .`

If you would like to use an existing installation of [MFEM](https://github.com/mfem/mfem/)
(outside of Spack), you can write your own host-config file by following the examples in the
host-configs directory and defining MFEM_DIR as well.

WARNING: The only MFEM build supported at the moment is the Makefile one (not the CMake one, yet).

License
-------

Serac is licensed under the BSD 3-Clause license,
(BSD-3-Clause or https://opensource.org/licenses/BSD-3-Clause).

Copyrights and patents in the Serac project are retained by contributors.
No copyright assignment is required to contribute to Serac.

See [LICENSE](https://github.com/LLNL/serac/blob/master/LICENSE),
[COPYRIGHT](https://github.com/LLNL/serac/blob/master/COPYRIGHT), and
[NOTICE](https://github.com/LLNL/serac/blob/master/NOTICE) for details.

Unlimited Open Source - BSD 3-clause Distribution
`LLNL-CODE-XXXXXX`  `OCEC-XX-XXX`

SPDX usage
------------

Individual files contain SPDX tags instead of the full license text.
This enables machine processing of license information based on the SPDX
License Identifiers that are available here: https://spdx.org/licenses/

Files that are licensed as BSD 3-Clause contain the following
text in the license header:

    SPDX-License-Identifier: (BSD-3-Clause)

External Packages
-----------------

Serac bundles some of its external dependencies in its repository.  These
packages are covered by various permissive licenses.  A summary listing
follows.  See the license included with each package for full details.


[//]: # (Note: The spaces at the end of each line below add line breaks)

PackageName: BLT  
PackageHomePage: https://github.com/LLNL/blt  
PackageLicenseDeclared: BSD-3-Clause  
