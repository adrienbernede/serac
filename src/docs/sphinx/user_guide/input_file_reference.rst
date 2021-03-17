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

Some keywords are not tied to a particular physics module - these are typically
global configuration options.

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


Solid Module Keywords
---------------------

The following keywords are available for configuring Serac's solid mechanics module:

=======================   ===========
Keyword                   Description
=======================   ===========
`boundary_conds`_         List of boundary conditions to apply.
`density`_                Material density.
`dynamics`_               Options for mass matrix inversion.
`equation_solver`_        Linear/nonlinear solver parameters for stiffness matrix.
`geometric_nonlin`_       Whether to include geometric nonlinearities.
`initial_displacement`_   The initial state of the displacement field.
`initial_velocity`_       The initial state of the velocity field.
`K`_                      Lamé's first parameter, bulk modulus.
`material_nonlin`_        Whether the material to model is linear elastic.
`mu`_                     Lamé's second parameter, shear modulus.
`order`_                  Order of the finite elements/polynomial interpolation.
`viscosity`_              Material viscosity.
=======================   ===========


Thermal Module Keywords
-----------------------

The following keywords are available for configuring Serac's thermal conduction module:

=======================   ===========
Keyword                   Description
=======================   ===========
`boundary_conds`_         List of boundary conditions to apply.
`cp`_                     Material specific heat capacity.
`dynamics`_               Options for mass matrix inversion.
`equation_solver`_        Linear/nonlinear solver parameters for stiffness matrix.
`initial_temperature`_    The initial state of the temperature field.
`initial_velocity`_       The initial state of the velocity field.
`kappa`_                  Material thermal conductivity.
`nonlinear_reaction`_     Options for a nonlinear reaction term.
`order`_                  Order of the finite elements/polynomial interpolation.
`rho`_                    Material density.
=======================   ===========


boundary_conds
--------------

Boundary conditions are defined as a dictionary - the key for a boundary condition must include the name
of the field on which it should be applied.  Note that only one of ``scalar_function``, ``vector_function``,
``constant``, ``vector_constant``, and ``piecewise_constant`` can be defined for a given boundary condition.
The following fields are available for each element of the dictionary:

**boundary_cond.attrs**

The mesh attributes on which the boundary condition should be applied.

.. code-block:: Lua

  thermal_conduction = {
    boundary_conds {
      ['temperature'] = {
        attrs = {1, 2, 3}
      }
    }
  }

**boundary_cond.component**

The vector component on which to apply a scalar coefficient.

.. code-block:: Lua

  thermal_conduction = {
    boundary_conds {
      ['temperature'] = {
        component = 1
      }
    }
  }

**boundary_cond.constant**

The scalar used to define an ``mfem::ConstantCoefficient``.

.. code-block:: Lua

  thermal_conduction = {
    boundary_conds {
      ['temperature'] = {
        constant = 2.5
      }
    }
  }

**boundary_cond.piecewise_constant**

The scalar used to define an ``mfem::PWConstCoefficient`` that maps mesh attributes to constants.

.. code-block:: Lua

  thermal_conduction = {
    boundary_conds {
      ['temperature'] = {
        piecewise_constant = {
          [1] = 3.0,
          [4] = 6.1
        }
      }
    }
  }

**boundary_cond.scalar_function**

The function used to define an ``mfem::FunctionCoefficient``.  Note that if the function is not
time-dependent, accepting a second parameter for time is not required.

.. code-block:: Lua

  thermal_conduction = {
    boundary_conds {
      ['temperature'] = {
        scalar_function = function(v, t)
          return v:norm() * t
        end
      }
    }
  }

**boundary_cond.vector_constant**

The scalar used to define an ``mfem::VectorConstantCoefficient``.

.. code-block:: Lua

  thermal_conduction = {
    boundary_conds {
      ['temperature'] = {
        vector_constant = {
          x = 0.0,
          y = 0.0,
          z = 0.0
        }
      }
    }
  }

**boundary_cond.vector_function**

The function used to define an ``mfem::VectorFunctionCoefficient``.  Note that if the function is not
time-dependent, accepting a second parameter for time is not required.

.. code-block:: Lua

  thermal_conduction = {
    boundary_conds {
      ['temperature'] = {
        vector_function = function(v, t)
          return v * t
        end
      }
    }
  }

**boundary_cond.vector_piecewise_constant**

The scalar used to define an ``mfem::VectorArrayCoefficient`` that maps mesh attributes to constants.

.. code-block:: Lua

  thermal_conduction = {
    boundary_conds {
      ['temperature'] = {
        vector_piecewise_constant = {
          [1] = {
            x = 1.0,
            y = 0.5,
            z = 0.0
          },
          [4] = {
            x = 2.5,
            y = 0.25,
            z = 1.0
          }
        }
      }
    }
  }


equation_solver
---------------

Wraps an ``mfem::Solver`` - underlying solver can be linear or nonlinear dependending on whether the ``nonlinear`` options are specified.

**equation_solver.linear.type**

Can be ``iterative`` or ``direct`` (for SuperLU).  ``linear.iterative_options`` is used when ``iterative`` is selected,
and ``linear.direct_options`` is used when ``direct`` is selected.

.. code-block:: Lua

  thermal_conduction = {
    equation_solver = {
      linear = {
        type = 'iterative'
      }
    }
  }
