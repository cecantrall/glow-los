"""N2 Vegard-Kaplan emission spectrum anchored to native GLOW output.

GLOW carries a single VK(0,9) volume-emission component rather than a
vibrationally resolved A-state model.  This module preserves that native
component exactly and distributes the other A-X bands using the Gilmore et al.
(1992) wavelengths and Einstein coefficients.

The default upper-state populations are a normalized geometric distribution.
Its first population is constrained by GLOW's own effective VK(0,9) rate,
``B(40) * A(10) = 0.0154 s^-1``, and the published A(0,9) coefficient.  This is
an explicit moderate-fidelity approximation, not a replacement for a cascade
and altitude-dependent quenching calculation.
"""

from __future__ import annotations

import csv
from functools import lru_cache
from importlib.resources import files

import numpy as np


GLOW_VK_09_COMPONENT = (1, 0)
GLOW_VK_09_EFFECTIVE_RATE_S = 0.02 * 0.7700
VK_09_UPPER = 0
VK_09_LOWER = 9
MAX_UPPER_VIBRATIONAL_LEVEL = 10


@lru_cache(maxsize=1)
def bands():
    """Return the bundled Gilmore et al. A-X band data as a structured array."""

    resource = files(__package__).joinpath("data/n2_vk_bands.csv")
    rows = []
    with resource.open("r", encoding="ascii", newline="") as stream:
        for row in csv.DictReader(stream):
            rows.append(
                (
                    int(row["v_upper"]),
                    int(row["v_lower"]),
                    float(row["wavelength_angstrom"]),
                    float(row["einstein_a_s-1"]),
                )
            )
    dtype = [
        ("v_upper", "i4"),
        ("v_lower", "i4"),
        ("wavelength_angstrom", "f8"),
        ("einstein_a_s", "f8"),
    ]
    return np.asarray(rows, dtype=dtype)


@lru_cache(maxsize=1)
def default_vibrational_populations():
    """Return the normalized upper-state populations used by default."""

    table = bands()
    anchor = table[
        (table["v_upper"] == VK_09_UPPER)
        & (table["v_lower"] == VK_09_LOWER)
    ]
    if anchor.size != 1:
        raise RuntimeError("The VK band table must contain one (0,9) anchor.")

    population_zero = GLOW_VK_09_EFFECTIVE_RATE_S / anchor["einstein_a_s"][0]
    ratio = _geometric_ratio(population_zero, MAX_UPPER_VIBRATIONAL_LEVEL + 1)
    populations = population_zero * ratio ** np.arange(
        MAX_UPPER_VIBRATIONAL_LEVEL + 1
    )
    populations /= populations.sum()
    populations.setflags(write=False)
    return populations


def from_glow(ceta, wavelengths, *, vibrational_populations=None):
    """Return a VK volume-emission spectrum derived from GLOW.

    Parameters
    ----------
    ceta
        Native GLOW component array with shape ``(process, feature, altitude)``.
    wavelengths
        Output bin centers in Angstroms.
    vibrational_populations
        Optional populations for A-state levels v'=0..10. Values are
        normalized internally. The relative population at v'=0 must be
        positive because native VK(0,9) is the absolute anchor.

    Returns
    -------
    numpy.ndarray
        Array with shape ``(altitude, wavelength)`` in ``cm^-3 s^-1`` per bin.
    """

    ceta = np.asarray(ceta, dtype=float)
    if ceta.ndim != 3 or ceta.shape[0] < 2 or ceta.shape[1] < 1:
        raise ValueError("ceta must contain the native GLOW VK(0,9) component.")
    wave = np.asarray(wavelengths, dtype=float)
    if wave.ndim != 1 or wave.size == 0 or not np.all(np.isfinite(wave)):
        raise ValueError("wavelengths must be a nonempty finite one-dimensional grid.")

    populations = _populations(vibrational_populations)
    table = bands()
    anchor_mask = (
        (table["v_upper"] == VK_09_UPPER)
        & (table["v_lower"] == VK_09_LOWER)
    )
    anchor_rate = table["einstein_a_s"][anchor_mask][0]
    anchor_profile = ceta[GLOW_VK_09_COMPONENT]
    spectrum = np.zeros((anchor_profile.size, wave.size), dtype=float)

    for band in table:
        if band["wavelength_angstrom"] < wave.min() or band["wavelength_angstrom"] > wave.max():
            continue
        relative_rate = (
            populations[band["v_upper"]]
            / populations[VK_09_UPPER]
            * band["einstein_a_s"]
            / anchor_rate
        )
        index = int(np.argmin(np.abs(wave - band["wavelength_angstrom"])))
        spectrum[:, index] += anchor_profile * relative_rate
    return spectrum


def _populations(values):
    if values is None:
        return default_vibrational_populations()
    populations = np.asarray(values, dtype=float)
    expected = MAX_UPPER_VIBRATIONAL_LEVEL + 1
    if populations.shape != (expected,):
        raise ValueError(f"vibrational_populations must contain {expected} values.")
    if not np.all(np.isfinite(populations)) or np.any(populations < 0.0):
        raise ValueError("vibrational_populations must be finite and nonnegative.")
    total = populations.sum()
    if total <= 0.0 or populations[0] <= 0.0:
        raise ValueError("vibrational_populations and its v'=0 value must be positive.")
    return populations / total


def _geometric_ratio(first_probability, count):
    if not 0.0 < first_probability < 1.0:
        raise ValueError("first_probability must be between zero and one.")
    lo, hi = 0.0, 1.0
    for _ in range(100):
        ratio = (lo + hi) / 2.0
        series_sum = (1.0 - ratio**count) / (1.0 - ratio)
        if 1.0 / series_sum > first_probability:
            lo = ratio
        else:
            hi = ratio
    return (lo + hi) / 2.0
