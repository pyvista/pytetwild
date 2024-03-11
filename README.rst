pytetwild
#########

|pypi| |MPL|

.. |pypi| image:: https://img.shields.io/pypi/v/pytetwild.svg?logo=python&logoColor=white
   :target: https://pypi.org/project/pytetwild/

.. |MPL| image:: https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg
   :target: https://opensource.org/license/mpl-2-0

``pytetwild`` is a Python library for mesh tetrahedralization. It is a
Python wrapper around the efficient C++ library for tetrahedral meshing provided by
`fTetWild <https://github.com/wildmeshing/fTetWild>`_.


Installation
************

We have pre-built wheels for Python 3.8 - Python 3.12 for Windows and Linux x64.

The recommended way to install ``pytetwild`` is via PyPI:

.. code:: sh

   pip install pytetwild

You can also clone the repository and install it from source, but since there's
C++ involved, the build is a bit more complicated. See ``CONTRIBUTING.md`` for
more details.


Usage
*****

To tetrahedralize a surface mesh from `PyVista <https://docs.pyvista.org>`_:

.. code:: py

   import pyvista as pv
   import pytetwild

   # Load or create a PyVista PolyData surface mesh
   # Here, we'll create a simple sphere mesh as an example
   surface_mesh = pv.Icosphere(nsub=2)

   # Convert the surface mesh to a tetrahedral mesh. For this example let's
   # use a coarse mesh
   tetrahedral_mesh = pytetwild.tetrahedralize_pv(surface_mesh, edge_length_fac=1))

   # Visualize the tetrahedral mesh in an "exploded" view
   tetrahedral_mesh.explode(1).plot(show_edges=True)

.. image:: https://github.com/pyvista/pytetwild/raw/main/exploded-sphere.png

You can also work with raw arrays. Here's a simple cube that we turn into tetrahedra.

.. code:: pycon

   import numpy as np

   # Define vertices of the cube
   vertices = np.array([
       [0, 0, 0],  # Vertex 0
       [1, 0, 0],  # Vertex 1
       [1, 1, 0],  # Vertex 2
       [0, 1, 0],  # Vertex 3
       [0, 0, 1],  # Vertex 4
       [1, 0, 1],  # Vertex 5
       [1, 1, 1],  # Vertex 6
       [0, 1, 1]   # Vertex 7
   ])

   # Define faces using vertex indices
   # Each face is a rectangle (also accepts triangles)
   faces = np.array([
       [0, 1, 2, 3],  # Front face
       [1, 5, 6, 2],  # Right face
       [5, 4, 7, 6],  # Back face
       [4, 0, 3, 7],  # Left face
       [4, 5, 1, 0],  # Bottom face
       [3, 2, 6, 7]   # Top face
   ])
   v_out, tetra = pytetwild.tetrahedralize(vertices, faces, optimize=False)


Usage - Options
---------------
We've surfaced a handful of parameters to each of our interfaces
``tetrahedralize`` and ``tetrahedralize_pv``. Here are the optional parameters.

.. code::

    Additional Parameters
    ---------------------
    edge_length_fac : float, default: 0.05
        Tetrahedral edge length as a function of bounding box diagional. The
        default ideal edge length is bb/20 (bounding box divided by 20).
    optimize : bool
        Improve the minimum scaled Jacobean for each cell. This leads to higher
        cell quality at the expense of computation time.


License and Acknowledgments
***************************

This project relies on ``fTetWild`` and credits go to the original authors for
their efficient C++ library for tetrahedral meshing. That work is licensed
under the Mozilla Public License v2.0.

The work in this repository is also licensed under the Mozilla Public License v2.0.

Support
*******

If you are having issues, please feel free to raise an `Issue
<https://github.com/pyvista/pytetwild/issues>`_.
