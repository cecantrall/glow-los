"""Approximate atomic-nitrogen 174.3 nm emission from GLOW output.

GLOW resolves the photodissociation and electron-impact contributions to the
N I 149.3 nm emission separately.  This module scales those two native
components to 174.3 nm.  The default ratio is the 200 eV cross-section ratio
reported by Mumma and Zipf (1971), 0.43 / 1.20.  Callers may supply different
ratios when a more suitable laboratory data set is available.
"""

from __future__ import annotations

import numpy as np


NI_1743_WAVELENGTH_ANGSTROM = 1743.0
NI_1493_PHOTODISSOCIATION_COMPONENT = (0, 13)
NI_1493_ELECTRON_IMPACT_COMPONENT = (1, 13)
NI_1743_TO_1493_RATIO = 0.43 / 1.20


def components_from_glow(
    ceta,
    *,
    photodissociation_ratio=NI_1743_TO_1493_RATIO,
    electron_impact_ratio=NI_1743_TO_1493_RATIO,
):
    """Return N I 174.3 nm production components in ``cm^-3 s^-1``.

    Parameters
    ----------
    ceta
        Native GLOW component array with shape ``(process, feature, altitude)``.
    photodissociation_ratio, electron_impact_ratio
        Ratios of N I 174.3 nm to 149.3 nm production for the corresponding
        native GLOW process.
    """

    ceta = np.asarray(ceta, dtype=float)
    if ceta.ndim != 3 or ceta.shape[0] < 2 or ceta.shape[1] <= 13:
        raise ValueError("ceta must contain the native GLOW N I 149.3 components.")

    photo_ratio = _nonnegative_scalar(
        photodissociation_ratio, "photodissociation_ratio"
    )
    electron_ratio = _nonnegative_scalar(
        electron_impact_ratio, "electron_impact_ratio"
    )
    photo = ceta[NI_1493_PHOTODISSOCIATION_COMPONENT] * photo_ratio
    electron = ceta[NI_1493_ELECTRON_IMPACT_COMPONENT] * electron_ratio
    return {"photodissociation": photo, "electron_impact": electron}


def from_glow(
    ceta,
    wavelengths,
    *,
    photodissociation_ratio=NI_1743_TO_1493_RATIO,
    electron_impact_ratio=NI_1743_TO_1493_RATIO,
):
    """Return the N I 174.3 nm volume-emission spectrum.

    The output has shape ``(altitude, wavelength)`` and units
    ``cm^-3 s^-1`` per wavelength bin.  The unresolved line is assigned to the
    nearest wavelength bin; an instrumental line shape can be applied later.
    """

    wave = _wavelength_grid(wavelengths)
    components = components_from_glow(
        ceta,
        photodissociation_ratio=photodissociation_ratio,
        electron_impact_ratio=electron_impact_ratio,
    )
    profile = components["photodissociation"] + components["electron_impact"]
    spectrum = np.zeros((profile.size, wave.size), dtype=float)
    index = int(np.argmin(np.abs(wave - NI_1743_WAVELENGTH_ANGSTROM)))
    spectrum[:, index] = profile
    return spectrum


def _nonnegative_scalar(value, name):
    value = float(value)
    if not np.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative.")
    return value


def _wavelength_grid(wavelengths):
    wave = np.asarray(wavelengths, dtype=float)
    if wave.ndim != 1 or wave.size == 0 or not np.all(np.isfinite(wave)):
        raise ValueError("wavelengths must be a nonempty finite one-dimensional grid.")
    return wave
