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

To tetrahedralize a surface mesh:

.. code:: pycon

   >>> import pytetwild
   >>> mesh = pytetwild.tetrahedralize("input_mesh.ply")
   >>> mesh.vertices
   array([[x1, y1, z1],
          [x2, y2, z2],
          ...], dtype=float32)
   >>> mesh.tetrahedra
   array([[i1, j1, k1, l1],
          [i2, j2, k2, l2],
          ...], dtype=int32)

You can also load a mesh, perform tetrahedralization, and export the tetrahedral mesh:

.. code:: pycon

    >>> import pytetwild
    >>> pytetwild.tetrahedralize("input_mesh.ply", "output_mesh.ply")



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
