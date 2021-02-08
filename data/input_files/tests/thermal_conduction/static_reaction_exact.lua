-- Comparison information
expected_t_l2norm = 0.643674
epsilon = 0.00001

-- Simulation time parameters
dt      = 1.0

main_mesh = {
    type = "file",
    -- mesh file
    mesh = "../../../meshes/star.mesh",
    -- serial and parallel refinement levels
    ser_ref_levels = 1,
    par_ref_levels = 1,
}

-- Solver parameters
thermal_conduction = {
    equation_solver = {
        linear = {
            type = "direct",
            direct_options = {
                print_level = 0,
            },
        },

        nonlinear = {
            rel_tol     = 1.0e-6,
            abs_tol     = 1.0e-12,
            max_iter    = 500,
            print_level = 1,
        },
    },

    -- polynomial interpolation order
    order = 3,

    -- material parameters
    kappa = 1.0,

    -- add a nonlinear reaction
    nonlinear_reaction = {
        reaction_function = function (temp)
            return temp^2 + 5.0
        end,
        d_reaction_function = function (temp)
            return 2.0 * temp
        end
    },

    source = {
        coef = function (v)
            return v.x^4 * v.y^2 - 2.0 * v.y + 5.0
        end
    },

    -- boundary condition parameters
    boundary_conds = {
        ['temperature'] = {
            attrs = {1},
            coef = function (v)
                return v.x^2 * v.y
            end
        },
    },
}