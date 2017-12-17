import numpy as np



def build_grid(geometric_properties, physical_properties, initial_condition, boundary_condition):
    # Aliases
    size_x = geometric_properties.size_x
    n_x = geometric_properties.n_x
    size_y = geometric_properties.size_y
    n_y = geometric_properties.n_y

    # Computed variables
    dx = size_x / n_x
    dy = size_y / n_y

    # Data among the grid
    G = {}
    G['dx'] = dx
    G['dy'] = dy
    G['n_x'] = n_x
    G['n_y'] = n_y
    G['k'] = physical_properties.k
    G['rho'] = physical_properties.rho
    G['c_p'] = physical_properties.c_p
    G['T'] = initial_condition.T

    return G
