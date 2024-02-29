"""Setup for pytetwild."""
import sys

from setuptools import find_packages

try:
    from skbuild import setup
except ImportError:
    print("Please update pip to pip 10 or greater, or a manually install the PEP 518 requirements in pyproject.toml", file=sys.stderr)
    raise

cmake_args = []
debug = False
cfg = 'Debug' if debug else 'Release'

setup(
    name='pytetwild',
    version='0.1.dev0',
    author='Alex Kaszynski',
    author_email='akascap@gmail.com',
    description='Python wrapper of fTetWild',
    long_description=open("README.rst").read(),
    long_description_content_type="text/rst",
    packages=find_packages('src'),
    package_dir={'':'src'},
    zip_safe=False,
    include_package_data=False,
    cmake_args=cmake_args=cmake_args + ['-DCMAKE_BUILD_TYPE=' + cfg],
    cmake_install_dir="src/",
    cmake_install_target='install',
    install_requires="numpy",
)
