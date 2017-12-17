from copy import deepcopy
import numpy as np
from petsc4py import PETSc
from solver.grid import build_grid
from solver.residual_function import residual_function


def create_solver(use_multigrid, use_newton, n_x, n_y):
    snes = PETSc.SNES().create()

    options = PETSc.Options()
    options.clear()
    if use_multigrid:
        options.setValue('pc_type', 'gamg')
        options.setValue('pc_gamg_reuse_interpolation', 'True')
    else:
        options.setValue('pc_type', 'ilu')
    options.setValue('ksp_type', 'gmres')
    if use_newton:
        options.setValue('snes_type', 'newtonls')
    else:
        options.setValue('snes_type', 'ksponly')
    options.setValue('ksp_atol', 1e-7)
    snes.setFromOptions()

    dmda = PETSc.DMDA().create([n_x, n_y], dof=1, stencil_width=1, stencil_type='star')
    snes.setDM(dmda)
    return snes


def solve(
    geometric_properties,
    physical_properties,
    initial_condition,
    boundary_condition,
    timestep_properties,

    # Output options
    on_timestep_callback=None,

    # Solver options
    use_multigrid=True,
    use_newton=True,
):
    if on_timestep_callback is None:
        on_timestep_callback = lambda *args, **kwargs: None

    grid = build_grid(
        geometric_properties,
        physical_properties,
        initial_condition,
        boundary_condition,
    )

    n_x = geometric_properties.n_x
    n_y = geometric_properties.n_y
    snes = create_solver(use_multigrid, use_newton, n_x, n_y)
    r = PETSc.Vec().createSeq(n_x * n_y)  # residual vector
    x = PETSc.Vec().createSeq(n_x * n_y)  # solution vector
    b = PETSc.Vec().createSeq(n_x * n_y)  # right-hand side

    current_time = 0.0
    old_grid = deepcopy(grid)
    on_timestep_callback(current_time, grid['T'])

    # solution = initial_condition
    while current_time < timestep_properties.final_time:
        snes.setFunction(residual_function, r, [current_time, grid, old_grid, timestep_properties, boundary_condition])
        initial_guess = np.array(old_grid['T'])
        x.setArray(initial_guess)
        b.set(0)
        snes.solve(b, x)
        solution = x[:].reshape(n_y, n_x)
        solution.flags['WRITEABLE'] = False

        # Retrieve the solution
        old_grid = deepcopy(grid)
        grid['T'][:] = solution

        # Advance in time
        current_time += timestep_properties.delta_t

        on_timestep_callback(current_time, solution)

    return grid
