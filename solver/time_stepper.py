from copy import deepcopy
import numpy as np
from solver.grid import build_grid
from solver.residual_function import residual_function
from collections import namedtuple


NonlinearSolverContext = namedtuple(
    'NonlinearSolverContext',
    [
        'current_time',
        'grid',
        'old_grid',
        'timestep_properties',
    ],
)

def transient_solve(
    geometric_properties,
    physical_properties,
    initial_condition,
    boundary_condition,
    timestep_properties,
    nonlinear_solver,    
    # Callbacks
    setup_nonlinear_solver=None,
    on_timestep_callback=None,
):
    if setup_nonlinear_solver is None:
        setup_nonlinear_solver = lambda nonlinear_solver, residual_function, context: None
    if on_timestep_callback is None:
        on_timestep_callback = lambda *args, **kwargs: None

    grid = build_grid(
        geometric_properties,
        physical_properties,
        initial_condition,
        boundary_condition,
    )

    current_time = 0.0
    old_grid = deepcopy(grid)
    on_timestep_callback(current_time, grid.T)

    while current_time < timestep_properties.final_time:
        context = NonlinearSolverContext(
            current_time,
            grid,
            old_grid,
            timestep_properties,
        )
        solve = setup_nonlinear_solver(
            nonlinear_solver,
            residual_function,
            context
        )

        initial_guess = np.array(old_grid.T)
        solution = solve(initial_guess)
        solution.flags['WRITEABLE'] = False

        # Retrieve the solution
        old_grid = deepcopy(grid)
        grid.T[:] = solution

        # Advance in time
        current_time += timestep_properties.delta_t
        on_timestep_callback(current_time, solution)

    return grid
