# Contributing to Our Project

Thank you for considering contributions to our project! Here's how you can help:

## Setting Up for Development

We're using [scikit-build-core](https://github.com/scikit-build/scikit-build-core) for our build system and building using cmake. You'll need a handful of system dependencies to build locally. On Linux, install them with:

```
sudo apt-get update
sudo apt-get install libblas-dev libboost-filesystem-dev libboost-system-dev libboost-thread-dev libglu1-mesa-dev libsuitesparse-dev xorg-dev ccache
```

To install the library in editable mode:

  1. Clone the repository:
     ```
     git clone https://github.com/pyvista/pytetwild
     ```
  2. Initialize submodules:
     ```
     git submodule update --init --recursive
     ```
  3. Install the project in editable mode and include development dependencies:
     ```
     pip install -e .[dev]
     ```

## Code Style and Quality

- **Documentation**: Follow the `numpydoc` style for docstrings.
- **Linting**: We use `ruff` for code styling to ensure consistency.
- **Pre-commit**: Utilize `pre-commit` hooks to automate checks before commits. Set up your local environment with:
  ```
  pip install pre-commit
  pre-commit install
  ```
- **Testing**: Write tests for new features or bug fixes and run them using `pytest`. Tests are located in the `tests` directory.

## How to Contribute

- **Bugs and Feature Requests**: Submit them through our [issues page](https://github.com/pyvista/pytetwild/issues).
- **Code Contributions**: Make changes in a fork of this repository and submit a pull request (PR) through our [PR page](https://github.com/pyvista/pytetwild/pulls). Follow the standard fork-and-PR workflow for contributions.

## Community and Conduct

- We expect all contributors to follow our [Code of Conduct](https://github.com/pyvista/pyvista/blob/main/CODE_OF_CONDUCT.md) to maintain a welcoming and inclusive community.

## License

- Contributions are subject to the project's license as found in the [LICENSE.md](./LICENSE.md) file.

We welcome your contributions and look forward to collaborating with you!
