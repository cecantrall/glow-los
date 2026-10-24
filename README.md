# glow-los

Calculate airglow radiances for **limb** and **disk** observations using GLOW,
REDISTER, and line-of-sight integration.

| Calculate | Use |
| --- | --- |
| O I 135.6 nm, N I 149.3 nm, and N2 LBH spectra | `glow_los.los` |
| Limb O1356, LBH short, and LBH long band radiances | `glow_los.limb_colors` |
| MSIS/IRI background atmospheres | `glow_los.background` |
| Spectral convolution and custom wavelength bands | `glow_los.spectral` |
| Optional Vegard–Kaplan and N I 174.3 nm emissions | [Extensions](docs/extensions.md) |

## Clone, install, build

```bash
git clone https://github.com/cecantrall/glow-los.git
cd glow-los
conda env create -f environment.yml
conda activate glow-los
python -m pip install -e '.[examples,dev]'
glow-los-build
```

The environment supplies Python 3.12, GNU Fortran, and NetCDF-Fortran.
Linux and macOS are the initial targets; on Windows use WSL.
Inside `caro-suite`, first run `cd packages/glow-los` and use the same commands.
No other `caro-suite` package is required.

Already have an environment? Install with `pip`, then follow the
[native build notes](docs/build.md).

## Run examples

```bash
python examples/glow_los_disk_spectra.py --oza 45 --no-show
python examples/glow_los_limb_spectra.py --tan-km 150 --no-show
```

Both examples calculate O1356 + N1493 + LBH at a fixed daytime location.
They print line radiances in Rayleighs and save a spectrum PNG and CSV in
`examples/output/`. OZA is measured from local zenith; tangent altitude is in km.
Remove `--no-show` to open the plot.

Change the instrument width or output folder:

```bash
python examples/glow_los_limb_spectra.py --tan-km 200 --psf-fwhm-a 3 --output-dir ./output --no-show
```

Edit `time` (UTC), `lat`, `lon`, `f107`, `f107a`, and `ap` in the scripts to
change the atmosphere. Wavelengths and PSF widths use Angstroms (10 Å = 1 nm).
The plotted spectrum is radiance per wavelength bin, not per Angstrom.

## Use in Python

After creating `bks`, `sza`, and `model` as in either example:

```python
from glow_los.spectral import to_rayleighs

limb, _, _ = model.los_1d("o1356", lat, lon, sza, bks, tan=150)
disk, _, _ = model.los_1d("o1356", lat, lon, sza, bks, oza=45)
print(to_rayleighs(limb), to_rayleighs(disk))
```

Specify one geometry per call. `los_1d` returns radiance, source, and path
contribution; LBH also returns its wavelength grid. Convert native radiances
with `to_rayleighs`; spectral combination helpers already return Rayleighs.
Use at most 99 altitude levels for REDISTER. The examples use 84.

## Checks and troubleshooting

```bash
python -m pytest tests -q
```

- Missing `glowpy` or `rt_fuv`: run `glow-los-build` with the active Python.
- Missing `netcdf.mod` or `libnetcdff`: see [build notes](docs/build.md).
- `IEEE_UNDERFLOW_FLAG`: rebuild; see [REDISTER notes](docs/redister.md).

GLOW is by Stan Solomon and REDISTER by Randy Gladstone. Preserve their source
notices; see the [GLOW](docs/glow.md) and [REDISTER](docs/redister.md) notes.
