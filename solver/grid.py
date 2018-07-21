from collections import namedtuple

Grid = namedtuple(
    'Grid',
    [
        'dx',
        'dy',
        'n_x',
        'n_y',
        'k',
        'rho',
        'c_p',
        'T',
        'boundary_condition'
    ],
)

def build_grid(geometric_properties, physical_properties, initial_condition, boundary_condition):
    # Aliases
    size_x = geometric_properties.size_x
    n_x = geometric_properties.n_x
    size_y = geometric_properties.size_y
    n_y = geometric_properties.n_y

    # Computed variables
    dx = size_x / n_x
    dy = size_y / n_y

    return Grid(
        dx,
        dy,
        n_x,
        n_y,
        physical_properties.k,
        physical_properties.rho,
        physical_properties.c_p,
        initial_condition.T,
        boundary_condition,
    )
