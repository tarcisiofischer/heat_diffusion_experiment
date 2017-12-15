import numpy as np
from time import time
from solver.properties import GeometricProperties, PhysicalProperties,\
    ConstantInitialCondition, GhostNodeBoundaryCondition, TimestepProperties
from solver.nonlinear_solver import solve


def plot_animated_results(result_list):
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    N_FRAMES = len(result_list)
    x, y = result_list[0].shape
    x = np.arange(0, x)
    y = np.arange(0, y)

    f = plt.figure()
    img = plt.imshow(result_list[0], vmin=np.min(result_list), vmax=np.max(result_list))

    def animate(i):
        img.set_data(result_list[i])

    anim = animation.FuncAnimation(
        f,
        animate,
        frames=N_FRAMES,
        interval=60,
        repeat_delay=1000
    )

    plt.show()


class InMemoryResultsHandler():

    def __init__(self):
        self._s = time()
        self.result_list = []

    def __call__(self, time_, result):
        print("Saving time %2.5f..." % (time_,), end='')
        self.result_list.append(np.array(result))
        print(" Done. (%s)" % (time() - self._s,))
        self._s = time()


def test_solver():
    result_handler = InMemoryResultsHandler()

    solve(
        GeometricProperties(
            n_x=100,
            n_y=100,
            size_x=1.0,
            size_y=1.0,
        ),
        PhysicalProperties(
            k=1.0,
            rho=1.0,
            c_p=1.0,
        ),
        ConstantInitialCondition(
            T=0.0,
        ),
        GhostNodeBoundaryCondition(
            T_E=lambda t: np.cos(t * 10.),
            T_W=lambda t: np.sin(t * 10.),
            T_N=lambda t: np.cos(t * 2.),
            T_S=lambda t: np.sin(t * 20.),
        ),
        TimestepProperties(
            delta_t=0.01,
            final_time=2.0,
        ),
        on_timestep_callback=result_handler
    )

    plot_animated_results(result_handler.result_list)
