import argparse
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from glow_los.background import bkmodel, solzen
from glow_los.los import los as LOSModel
from glow_los.spectral import convolved_emission_spectrum, to_rayleighs


def update_glow_params(model, f107=120, f107a=120, f107p=120, lbh_scale_fac=1.351):
    model.d_glow.update({"f107": f107, "f107a": f107a, "f107p": f107p, "lbh_scale_fac": lbh_scale_fac})


def parse_args():
    parser = argparse.ArgumentParser(description="Run a disk-viewing GLOW LOS spectrum example.")
    parser.add_argument("--oza", type=float, default=45.0, help="Observing zenith angle in degrees (0 <= OZA < 90).")
    parser.add_argument("--no-show", action="store_true", help="Save outputs without opening a plot window.")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "output")
    parser.add_argument(
        "--psf-fwhm-a",
        type=float,
        default=3.0,
        help="Gaussian spectral PSF full-width at half maximum in Angstroms.",
    )
    args = parser.parse_args()
    if not np.isfinite(args.oza) or not 0 <= args.oza < 90:
        parser.error("--oza must be between 0 (inclusive) and 90 (exclusive).")
    if not np.isfinite(args.psf_fwhm_a) or args.psf_fwhm_a <= 0:
        parser.error("--psf-fwhm-a must be positive and finite.")
    return args


def main():
    args = parse_args()
    time = datetime(2002, 12, 31, 0)
    lat = 0
    lon = 180
    ozas = [args.oza]
    f107 = 120
    f107a = 120
    ap = 1

    levs = range(80, 500, 5)
    bks = bkmodel.bks_from_msis_iri(time, lat, lon, levs, f107=f107, f107a=f107a, ap=ap)

    idate = solzen.datetime_to_idate(time)
    sza = solzen.solzen(idate, time.hour, lat, lon)

    model = LOSModel()
    update_glow_params(model, f107=f107, f107a=f107a)

    o1356 = []
    n1493 = []
    lbh = []
    wave = None
    for oza in ozas:
        toa1356, _, _ = model.los_1d("o1356", lat, lon, sza, bks, oza=oza)
        toalbh, _, _, this_wave = model.los_1d("n2lbh", lat, lon, sza, bks, oza=oza)
        toa1493, _, _ = model.los_1d("n1493", lat, lon, sza, bks, oza=oza)

        o1356.append(toa1356)
        n1493.append(toa1493)
        lbh.append(toalbh)
        wave = this_wave

    o1356 = np.asarray(o1356)
    n1493 = np.asarray(n1493)
    lbh = np.asarray(lbh)
    wave = np.asarray(wave)

    sensor_spectrum = convolved_emission_spectrum(
        wave,
        lbh[0],
        {1356.0: o1356[0], 1493.0: n1493[0]},
        args.psf_fwhm_a,
    )

    scale = 1e-3

    fig, ax = plt.subplots(1, 1, figsize=(8.5, 3))
    ax.plot(wave, sensor_spectrum * scale, "-", label=f"Combined spectrum, FWHM={args.psf_fwhm_a:g} A")
    ax.set_title("Disk Radiance Spectrum")
    ax.grid()
    ax.set_xlabel("Wavelength [A]")
    ax.set_ylabel("Radiance [kR / bin]")
    ax.legend()
    fig.tight_layout()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "glow_los_disk_spectra.png"
    fig.savefig(output_path, dpi=150)
    print(f"Wrote: {output_path}")
    csv_path = output_path.with_suffix(".csv")
    np.savetxt(csv_path, np.column_stack([wave, sensor_spectrum]), delimiter=",",
               header="wavelength_A,radiance_R_per_bin", comments="")
    print(f"Wrote: {csv_path}")
    print(
        "Disk radiances: "
        f"O1356={float(to_rayleighs(o1356[0])):.3f} R, "
        f"N1493={float(to_rayleighs(n1493[0])):.3f} R, "
        f"convolved spectrum max={float(np.nanmax(sensor_spectrum)):.3f} R"
    )

    if not args.no_show and matplotlib.get_backend().lower() != "agg":
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
