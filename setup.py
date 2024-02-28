"""Setup for pyminiply."""
from io import open as io_open
import os
import sys

from Cython.Build import cythonize
import numpy as np
from setuptools import Extension, setup

filepath = os.path.dirname(__file__)

# Define macros for cython
macros = []
if os.name == "nt":  # windows
    extra_compile_args = ["/O2", "/w", "/GS"]
elif os.name == "posix":  # linux org mac os
    if sys.platform == "linux":
        extra_compile_args = ["-std=gnu++11", "-O3", "-w"]
    else:  # probably mac os
        extra_compile_args = ["-std=c++11", "-O3", "-w"]
else:
    raise Exception(f"Unsupported OS {os.name}")


# Get version from version info
__version__ = None
version_file = os.path.join(filepath, "pyminiply", "_version.py")
with io_open(version_file, mode="r") as fd:
    exec(fd.read())

# readme file
readme_file = os.path.join(filepath, "README.rst")


setup(
    name="pytetwild",
    packages=["pytetwild"],
    version=__version__,
    description="Tetrahedralize surfaces using fTetWild",
    long_description=open(readme_file).read(),
    long_description_content_type="text/x-rst",
    author="PyVista Developers",
    author_email="info@pyvista.org",
    license="MPLv2",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.8",
    url="https://github.com/pyvista/pytetwild",
    # Build cython modules
    ext_modules=cythonize(
        [
            Extension(
                "pytetwild._wrapper",
                [
                    "pytetwild/src/FTetWildWrapper.cpp",
                    "pytetwild/_wrapper.pyx",
                    "pytetwild/fTetWild/src.cpp",
                ],
                language="c++",
                extra_compile_args=extra_compile_args,
                define_macros=macros,
                include_dirs=[np.get_include()],
            )
        ]
    ),
    package_data={
        "pytetwild": ["*.pyx", "*.hpp"],
        "pytetwild/wrapper": ["*.c", "*.h"],
    },
    keywords="tetrahedralize tetwild ftetwild",
    install_requires=["numpy>1.11.0"],
)
