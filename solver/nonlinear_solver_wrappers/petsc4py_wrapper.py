from petsc4py import PETSc
import functools

def create(use_multigrid=False, use_newton=True):
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

    return snes


def setup(nonlinear_solver, residual_function, context):
    dmda = PETSc.DMDA().create([context.grid.n_x, context.grid.n_y], dof=1, stencil_width=1, stencil_type='star')
    nonlinear_solver.setDM(dmda)

    n_x = context.grid.n_x
    n_y = context.grid.n_y
    r = PETSc.Vec().createSeq(n_x * n_y)  # residual vector

    def wrapped_residual_function(snes, X, f, t, grid, old_grid, timestep_properties):
        residual_function(
            snes,
            X.getArray(readonly=True).reshape((n_y, n_x)),
            f,
            t,
            grid,
            old_grid,
            timestep_properties
        )
    nonlinear_solver.setFunction(wrapped_residual_function, r, context)

    def wrapped_nonlinear_solver_callable(initial_guess):
        x = PETSc.Vec().createSeq(n_x * n_y)  # solution vector
        b = PETSc.Vec().createSeq(n_x * n_y)  # right-hand side
        x.setArray(initial_guess)
        b.set(0)
        nonlinear_solver.solve(b, x)
        solution = x[:].reshape(n_y, n_x)
        return solution

    return wrapped_nonlinear_solver_callable
