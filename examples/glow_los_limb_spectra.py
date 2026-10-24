import argparse
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from glow_los.background import bkmodel, solzen
from glow_los.los import los as LOSModel
from glow_los.spectral import convolved_emission_spectrum, to_rayleighs


def parse_args():
    parser = argparse.ArgumentParser(description="Run a limb-viewing GLOW LOS spectrum example.")
    parser.add_argument("--tan-km", type=float, default=150.0, help="Tangent altitude in km (80 <= altitude < 495).")
    parser.add_argument("--no-show", action="store_true", help="Save outputs without opening a plot window.")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "output")
    parser.add_argument(
        "--psf-fwhm-a",
        type=float,
        default=3.0,
        help="Gaussian spectral PSF full-width at half maximum in Angstroms.",
    )
    args = parser.parse_args()
    if not np.isfinite(args.tan_km) or not 80 <= args.tan_km < 495:
        parser.error("--tan-km must be within the example atmosphere: 80 <= altitude < 495 km.")
    if not np.isfinite(args.psf_fwhm_a) or args.psf_fwhm_a <= 0:
        parser.error("--psf-fwhm-a must be positive and finite.")
    return args


def main():
    args = parse_args()
    time = datetime(2002, 12, 31, 0)
    lat = 0
    lon = 180
    tanalts = [args.tan_km]
    f107 = 120
    f107a = 120
    ap = 1

    levs = range(80, 500, 5)
    bks = bkmodel.bks_from_msis_iri(time, lat, lon, levs, f107=f107, f107a=f107a, ap=ap)

    idate = solzen.datetime_to_idate(time)
    sza = solzen.solzen(idate, time.hour, lat, lon)

    model = LOSModel()
    model.d_glow.update({"f107": f107, "f107a": f107a})

    toa1356, _, _ = model.los_1d("o1356", lat, lon, sza, bks, tan=tanalts[0])
    toalbh, _, _, wave = model.los_1d("n2lbh", lat, lon, sza, bks, tan=tanalts[0])
    toa1493, _, _ = model.los_1d("n1493", lat, lon, sza, bks, tan=tanalts[0])

    sensor_spectrum = convolved_emission_spectrum(
        wave,
        toalbh,
        {1356.0: toa1356, 1493.0: toa1493},
        args.psf_fwhm_a,
    )

    scale = 1e-3
    fig, ax = plt.subplots(1, 1, figsize=(8.5, 3))
    ax.plot(wave, sensor_spectrum * scale, "-", label=f"Combined spectrum, FWHM={args.psf_fwhm_a:g} A")
    ax.set_title("Limb Radiance Spectrum")
    ax.grid()
    ax.set_xlabel("Wavelength [A]")
    ax.set_ylabel("Radiance [kR / bin]")
    ax.legend()
    fig.tight_layout()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "glow_los_limb_spectra.png"
    fig.savefig(output_path, dpi=150)
    print(f"Wrote: {output_path}")
    csv_path = output_path.with_suffix(".csv")
    np.savetxt(csv_path, np.column_stack([wave, sensor_spectrum]), delimiter=",",
               header="wavelength_A,radiance_R_per_bin", comments="")
    print(f"Wrote: {csv_path}")
    print(
        "Limb radiances: "
        f"O1356={float(to_rayleighs(toa1356)):.3f} R, "
        f"N1493={float(to_rayleighs(toa1493)):.3f} R, "
        f"convolved spectrum max={float(np.nanmax(sensor_spectrum)):.3f} R"
    )

    if not args.no_show and matplotlib.get_backend().lower() != "agg":
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
