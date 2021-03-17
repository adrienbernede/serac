.. ## Copyright (c) 2019-2021, Lawrence Livermore National Security, LLC and
.. ## other Serac Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

====================
Input File Reference
====================

This document contains a keyword reference for input files to Serac's driver.
The input files discussed in the :ref:`simple-conduction-label` tutorial use
a subset of these keywords as well.

Global Keywords
---------------

Some keywords are not tied to a particular physics module - typically these
are used for global configuration options or for the 

==============    ===========
Keyword           Description
==============    ===========
`dt`_             Simulation time step.
`main_mesh`_      The primary mesh used by all physics modules.
`output_type`_    The format to save simulation data in.
`t_final`_        The time to which the simulation should be run.
==============    ===========


dt
--

The time increment for each iteration of the simulation main loop.  Unused for quasistatic simulations.

.. code-block:: Lua

  dt = 0.05

main_mesh
---------

The problem mesh - can be read from a file, or generated automatically, in the case of regular geometrical
object meshes.

**main_mesh.approx_elements**

The approximate number of elements to generate, only applicable when **main_mesh.type** is ``disk`` or ``ball``.

.. code-block:: Lua

  main_mesh = {
    approx_elements = 1000
  }

**main_mesh.elements**

The number of elements to produce in each dimension, only applicable when **main_mesh.type** is ``box``.

.. code-block:: Lua

  main_mesh = {
    elements = {
      -- Produces a cuboid mesh with 20 x 30 x 10 = 6000 elements
      x = 20,
      y = 30,
      z = 10
    }
  }

**main_mesh.mesh**

The path to the mesh file, only applicable when **main_mesh.type** is ``file``.  If ``mesh`` is a relative
path, Serac's driver will search relative to the directory the input file is located in.

.. code-block:: Lua

  main_mesh = {
    mesh = '../meshes/beam-hex.mesh'
  }

**main_mesh.par_ref_levels**

The number of refinements to perform *after* distributing the mesh.

.. code-block:: Lua

  main_mesh = {
    par_ref_levels = 2
  }

**main_mesh.ser_ref_levels**

The number of refinements to perform *before* distributing the mesh.

.. code-block:: Lua

  main_mesh = {
    ser_ref_levels = 2
  }

**main_mesh.size**

The size of the mesh to generate, only applicable when **main_mesh.type** is ``box``.

.. code-block:: Lua

  main_mesh = {
    size = {
      -- Produces a cuboid mesh measuring 2 x 2 x 3
      x = 2,
      y = 2,
      z = 3
    }
  }

**main_mesh.type**

The type of the mesh to create.  For generated meshes, ``disk``, ``ball``, and ``box`` options
are available.  Use ``file`` to read a mesh from a file.

.. code-block:: Lua

  main_mesh = {
    type = 'disk'
  }


output_type
-----------

Serac supports the following output formats:

- `GLVis <https://glvis.org/>`_
- `ParaView <https://www.paraview.org/>`_
- `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit>`_
- `SidreVisIt <https://axom.readthedocs.io/en/develop/axom/sidre/docs/sphinx/mfem_sidre_datacollection.html>`_

.. code-block:: Lua

  output_type = 'GLVis'


t_final
-------

The point at which the simulation should be terminated.  Unused for quasistatic simulations.
Note that all simulations currently start at ``t = 0``, including restart runs.

.. code-block:: Lua

 t_final = 5
