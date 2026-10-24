from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from . import PACKAGE_ROOT


GLOW_SOURCES = (
    "cglow.f90 glowdriver.f90 glow.f90 bands.f90 egrid.f90 ephoto.f90 "
    "etrans.f90 exsect.f fieldm.f gchem.f90 maxt.f90 qback.f90 rcolum.f90 ssflux.f90"
).split()
RT_SOURCES = "rt_utils.f90 rt_lbh.f90 rt_o1356.f90 rt_n1493.f90".split()


def build() -> None:
    """Compile with the active Python and GNU Fortran toolchain."""
    env = os.environ.copy()
    compiler = env.get("FC", "gfortran")
    if shutil.which(compiler) is None:
        raise RuntimeError(f"Fortran compiler {compiler!r} was not found. Set FC to its executable path.")
    env["FC"] = compiler
    env["CC"] = env.get("CC", "cc")
    config = shutil.which("nf-config")
    for key, option, filename in (
        ("NETCDFF_INCLUDE", "--includedir", "netcdf.mod"),
        ("NETCDFF_LIB", "--prefix", None),
    ):
        if not env.get(key) and config:
            value = subprocess.check_output([config, option], text=True).strip()
            env[key] = str(Path(value) / "lib") if filename is None else value
        path = Path(env.get(key, "/nonexistent"))
        found = (path / filename).is_file() if filename else any(path.glob("libnetcdff.*"))
        if not found:
            raise RuntimeError(f"Set {key} to the directory containing {filename or 'libnetcdff'}.")

    # Meson is supported on Python 3.12+, where distutils is unavailable.
    env["FFLAGS"] = "-O2 -ffree-line-length-none -std=legacy " + env.get("FFLAGS", "")
    for directory, module, sources in (
        ("glow", "glowpy", GLOW_SOURCES), ("rt", "rt_fuv", RT_SOURCES),
    ):
        subprocess.run(
            [sys.executable, "-m", "numpy.f2py", "-c", "--backend", "meson",
             "-m", module, *sources, f"-I{env['NETCDFF_INCLUDE']}",
             f"-L{env['NETCDFF_LIB']}", "-lnetcdff"],
            cwd=PACKAGE_ROOT / directory, env=env, check=True,
        )
    # Underflow in optically thick paths is expected. Keep serious IEEE
    # diagnostics visible without changing floating-point arithmetic.
    subprocess.run(
        [compiler, "-g", "-fbacktrace", "-std=legacy",
         "-ffpe-summary=invalid,zero,overflow", *shlex.split(os.environ.get("FFLAGS", "")),
         "redister.f", "-o", "a.out"],
        cwd=PACKAGE_ROOT / "rt" / "redister_src", env=env, check=True,
    )


def main() -> None:
    try:
        build()
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        raise SystemExit(f"glow-los-build: {exc}") from exc
    print("Built GLOW, FUV transfer, and REDISTER. Run an example to check the installation.")


if __name__ == "__main__":
    main()
