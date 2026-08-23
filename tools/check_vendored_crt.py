"""Fail if a repaired wheel vendors a newer VC++ runtime than it was built with.

delvewheel resolves ``msvcp140.dll`` and ``concrt140.dll`` off ``PATH``. On
this project's Windows job a conda env is active, so they come from
conda-forge's ``vc14_runtime`` rather than from the toolset that compiled the
extension, and their version is then whatever conda-forge published most
recently -- an input nothing in this repo controls.

delvewheel does check the pair, but only warns when the redist it found is
*older* than the compiling toolset (see ``_toolset_too_old``). The conda case
is the other direction and passes silently. That is the gap this closes.

The mismatch is not cosmetic. ``vcruntime140.dll`` is *not* vendored -- CPython
ships it, so delvewheel treats it as a system DLL -- while ``msvcp140`` and
``concrt140`` are vendored under mangled names. A wheel that mixes toolsets is
therefore running one half of a redistributable set against the other, and the
mangled names mean the vendored copies load *alongside* the system ones rather
than instead of them, putting two Concurrency Runtimes in a single process.

Usage::

    python tools/check_vendored_crt.py <dir-or-wheel> [...]
"""

from __future__ import annotations

import re
import sys
import zipfile
from pathlib import Path

import pefile

# delvewheel appends "-<32 hex digits>" when it mangles a vendored DLL name.
MANGLE_SUFFIX = r"(-[0-9a-f]{32})?"

# The redistributable set, as delvewheel itself defines it in _dll_list.py.
VC_REDIST = re.compile(
    r"((vcruntime\d.*)|(vccorlib\d.*)|(msvcp[\d_].*)|(msvcr(t|\d.*))|(concrt\d.*)"
    r"|(mfc(m(ifc)?)?\d.*)|(vcamp\d.*)|(vcomp(\d.*|))|ucrtbase(d|_.*|))" + MANGLE_SUFFIX + r"\.dll"
)

EXTENSIONS = (".pyd", ".dll")


def linker_version(data: bytes) -> tuple[int, int]:
    """Return the ``(major, minor)`` linker version stamped in a PE file."""
    pe = pefile.PE(data=data, fast_load=True)
    try:
        return (
            pe.OPTIONAL_HEADER.MajorLinkerVersion,
            pe.OPTIONAL_HEADER.MinorLinkerVersion,
        )
    finally:
        pe.close()


def _bucket(version: tuple[int, int]) -> tuple[int, int]:
    """Collapse a linker version to the granularity a redist ships at.

    Matches delvewheel's own comparison: the minor version moves within a
    toolset without a redistributable actually changing, so compare
    ``minor // 10`` rather than ``minor``.
    """
    return version[0], version[1] // 10


def check_wheel(path: Path) -> list[str]:
    """Return a list of problems found in one wheel. Empty means it is clean."""
    built_by: dict[str, tuple[int, int]] = {}
    vendored: dict[str, tuple[int, int]] = {}

    with zipfile.ZipFile(path) as zf:
        for name in zf.namelist():
            if not name.lower().endswith(EXTENSIONS):
                continue
            base = name.rsplit("/", 1)[-1]
            try:
                version = linker_version(zf.read(name))
            except pefile.PEFormatError:
                continue
            if VC_REDIST.fullmatch(base.lower()):
                vendored[base] = version
            elif name.lower().endswith(".pyd"):
                built_by[base] = version

    if not built_by:
        return [f"{path.name}: no extension module found to compare against"]
    if not vendored:
        return []

    # The newest toolset that produced any of our own binaries.
    toolset = max(built_by.values())
    problems = [
        f"{path.name}: {base} was built by toolset {ver[0]}.{ver[1]}, newer than "
        f"the {toolset[0]}.{toolset[1]} that produced "
        f"{max(built_by, key=lambda k: built_by[k])}"
        for base, ver in sorted(vendored.items())
        if _bucket(ver) > _bucket(toolset)
    ]
    return problems


def main(argv: list[str]) -> int:
    """Check every wheel named on the command line."""
    if not argv:
        print(__doc__, file=sys.stderr)
        return 2

    wheels: list[Path] = []
    for arg in argv:
        target = Path(arg)
        wheels.extend(sorted(target.glob("*.whl")) if target.is_dir() else [target])

    if not wheels:
        print(f"No wheels found in {', '.join(argv)}", file=sys.stderr)
        return 1

    problems: list[str] = []
    for wheel in wheels:
        found = check_wheel(wheel)
        print(f"{'FAIL' if found else 'ok  '}  {wheel.name}")
        problems.extend(found)

    for problem in problems:
        print(f"  {problem}", file=sys.stderr)
    if problems:
        print(
            "\nThe vendored VC++ runtime is newer than the toolset that compiled "
            "the extension. delvewheel resolves it off PATH, so an active conda "
            "env supplies it instead of the compiler. Pin vc14_runtime to the "
            "toolset, or keep the conda env off PATH during the repair.",
            file=sys.stderr,
        )
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
