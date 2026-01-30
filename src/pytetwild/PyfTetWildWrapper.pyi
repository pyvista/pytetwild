import numpy as np
from numpy.typing import NDArray

def tetrahedralize_mesh(
    vertices: NDArray[np.float64],
    faces: NDArray[np.int32],
    optimize: bool,
    skip_simplify: bool,
    edge_length_r: float,
    epsilon: float,
    stop_energy: float,
    coarsen: bool,
    num_threads: int,
    loglevel: int,
    quiet: bool,
) -> tuple[np.ndarray, np.ndarray]: ...
def tetrahedralize_csg(
    csg_file: str,
    epsilon: float,
    edge_length_r: float,
    stop_energy: float,
    coarsen: bool,
    num_threads: int,
    log_level: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]: ...
