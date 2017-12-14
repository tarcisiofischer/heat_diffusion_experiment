import numpy as np
from solver.linear_solver import solve, GeometricProperties, PhysicalProperties, \
    ConstantInitialCondition, TimestepProperties, GhostNodeBoundaryCondition


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


def test_solver():
    result_list = []
    def save_results_in_memory(time, result):
        print("Saving time %2.5f..." % (time,), end='')
        result_list.append(np.array(result))
        print(" Done.")

    solve(
        GeometricProperties(
            n_x=20,
            n_y=20,
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
            T_E=3.0,
            T_W=3.0,
            T_N=3.0,
            T_S=3.0,
        ),
        TimestepProperties(
            delta_t=0.001,
            final_time=0.1,
        ),
        on_timestep_callback=save_results_in_memory
    )

    plot_animated_results(result_list)
