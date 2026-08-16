"""Create a requirements_testing.txt file directly from a pyproject.toml."""

import tomli

with open("pyproject.toml", "rb") as f:
    data = tomli.load(f)
deps = data["project"]["optional-dependencies"]["tests"]
with open("requirements-tests.txt", "w") as f:
    f.write("\n".join(deps))
