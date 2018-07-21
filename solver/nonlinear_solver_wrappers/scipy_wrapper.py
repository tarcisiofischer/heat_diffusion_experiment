import numpy as np
from scipy.optimize.nonlin import newton_krylov

def create():
    return newton_krylov


def setup(nonlinear_solver, residual_function, context):
    def residual_function_wrapped(x):
        f = np.zeros(shape=(context.grid.n_x * context.grid.n_y))
        residual_function(
            None,
            x,
            f,
            context.current_time,
            context.grid,
            context.old_grid,
            context.timestep_properties,
        )
        return f
    
    def wrapped_nonlinear_solver(initial_guess):
        return nonlinear_solver(residual_function_wrapped, initial_guess)
    
    return wrapped_nonlinear_solver
