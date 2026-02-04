import numpy as np
from numpy.typing import NDArray

def tetrahedralize_mesh(
    vertices: NDArray[np.float64],
    faces: NDArray[np.uint32],
    optimize: bool,
    skip_simplify: bool,
    edge_length_r: float,
    edge_length_abs: float,
    epsilon: float,
    stop_energy: float,
    coarsen: bool,
    num_threads: int,
    num_opt_iter: int,
    loglevel: int,
    quiet: bool,
    vtk_ordering: bool,
    disable_filtering: bool,
) -> tuple[NDArray[np.float64], NDArray[np.int32]]: ...
def tetrahedralize_csg(
    csg_file: str,
    epsilon: float,
    edge_length_r: float,
    stop_energy: float,
    coarsen: bool,
    num_threads: int,
    log_level: int,
    vtk_ordering: bool,
) -> tuple[NDArray[np.float64], NDArray[np.int32], NDArray[np.int32]]: ...
