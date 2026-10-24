# Native build

`glow-los-build` builds GLOW, FUV transfer, and the REDISTER executable using
the active Python. Run it again after changing Python, NumPy, or Fortran sources.

Required: GNU Fortran, a C compiler, NetCDF-Fortran, NumPy, Meson, and Ninja.
The recommended `environment.yml` supplies these.

The builder discovers NetCDF with `nf-config`. To select it explicitly:

```bash
export NETCDFF_INCLUDE="$CONDA_PREFIX/include"
export NETCDFF_LIB="$CONDA_PREFIX/lib"
glow-los-build
```

`NETCDFF_INCLUDE` must contain `netcdf.mod`; `NETCDFF_LIB` must contain
`libnetcdff`. Use a NetCDF library built with a compatible Fortran compiler.
Set `FC` to your GNU Fortran executable and `CC` to your C compiler if needed.
Do not copy `.so` files between machines or Python environments.

For a server without a display:

```bash
MPLBACKEND=Agg python examples/glow_los_disk_spectra.py --no-show
```

REDISTER writes temporary files under the system temporary directory, then
removes them. It does not need write access to another installed package.
