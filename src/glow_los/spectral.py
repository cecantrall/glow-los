from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from scipy.ndimage import gaussian_filter1d


def to_rayleighs(value):
    return np.asarray(value) / 1e6


def combine_emissions(wavelengths, lbh_spectrum, line_radiances: Mapping[float, float], line_fwhm=None):
    spectrum = np.asarray(lbh_spectrum, dtype=float).copy()
    wave = validate_wavelength_grid(wavelengths)
    validate_spectrum_grid(wave, spectrum)
    if line_fwhm is None:
        for line_wavelength, radiance in line_radiances.items():
            idx = int(np.argmin(np.abs(wave - float(line_wavelength))))
            spectrum[..., idx] += np.asarray(radiance, dtype=float)
    else:
        spectrum = spectrum + gaussian_line_spectrum(wave, line_radiances, line_fwhm)
    return to_rayleighs(spectrum)


def combine_contribution_spectra(wavelengths, lbh_contribution, line_contributions: Mapping[float, np.ndarray]):
    contribution = np.asarray(lbh_contribution, dtype=float).copy()
    wave = np.asarray(wavelengths, dtype=float)
    for line_wavelength, line_contribution in line_contributions.items():
        idx = int(np.argmin(np.abs(wave - float(line_wavelength))))
        contribution[:, idx] += np.asarray(line_contribution, dtype=float)
    return to_rayleighs(contribution)


def integrate_bands(wavelengths, spectrum, bands: Mapping[str, tuple[float, float]]):
    wave = validate_wavelength_grid(wavelengths)
    spectrum = np.asarray(spectrum)
    validate_spectrum_grid(wave, spectrum)
    values = {}
    for name, (lo, hi) in bands.items():
        mask = (wave >= lo) & (wave <= hi)
        values[name] = np.sum(spectrum[..., mask], axis=-1)
    return values


def convolve_spectrum(wavelengths, spectrum, fwhm, spacing=None):
    wave = validate_wavelength_grid(wavelengths)
    fwhm = _validate_positive_fwhm(fwhm)
    spectrum = np.asarray(spectrum, dtype=float)
    validate_spectrum_grid(wave, spectrum)
    if spacing is None:
        spacing = float(np.median(np.diff(wave)))
    spacing = float(spacing)
    if not np.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("wavelength spacing must be positive.")
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return gaussian_filter1d(spectrum, sigma / spacing, axis=-1, mode="nearest")


def gaussian_line_spectrum(wavelengths, line_radiances: Mapping[float, float], fwhm):
    """Return line radiances distributed onto `wavelengths` by a Gaussian PSF."""

    wave = validate_wavelength_grid(wavelengths)
    fwhm = _validate_positive_fwhm(fwhm)
    spectrum = None
    for line_wavelength, radiance in line_radiances.items():
        radiance = np.asarray(radiance, dtype=float)
        weights = gaussian_psf_weights(wave, line_wavelength, fwhm)
        if spectrum is None:
            spectrum = np.zeros(radiance.shape + (wave.size,), dtype=float)
        spectrum = spectrum + radiance[..., None] * weights
    if spectrum is None:
        return np.zeros(wave.shape, dtype=float)
    return spectrum


def gaussian_psf_weights(wavelengths, center, fwhm):
    """Return unit-sum Gaussian PSF weights on a wavelength grid."""

    wave = validate_wavelength_grid(wavelengths)
    fwhm = _validate_positive_fwhm(fwhm)
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    weights = np.exp(-0.5 * ((wave - float(center)) / sigma) ** 2)
    total = float(np.sum(weights))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("Gaussian PSF has zero support on the wavelength grid.")
    return weights / total


def convolved_emission_spectrum(wavelengths, lbh_spectrum, line_radiances: Mapping[float, float], psf_fwhm):
    """Return a combined LBH and line spectrum after Gaussian PSF convolution.

    Inputs use the native GLOW/RT radiance units. The returned spectrum is in
    Rayleighs on the input wavelength grid.
    """

    wave = validate_wavelength_grid(wavelengths)
    psf_fwhm = _validate_positive_fwhm(psf_fwhm)
    lbh_convolved = convolve_spectrum(wave, to_rayleighs(lbh_spectrum), psf_fwhm)
    line_convolved = to_rayleighs(gaussian_line_spectrum(wave, line_radiances, psf_fwhm))
    return lbh_convolved + line_convolved


def rebin_spectrum(wavelengths, spectrum, output_wavelengths):
    return np.interp(np.asarray(output_wavelengths, dtype=float), np.asarray(wavelengths, dtype=float), np.asarray(spectrum, dtype=float))


def convolve_and_rebin(wavelengths, spectrum, output_wavelengths, fwhm):
    return rebin_spectrum(wavelengths, convolve_spectrum(wavelengths, spectrum, fwhm), output_wavelengths)


def validate_wavelength_grid(wavelengths):
    """Return a finite, strictly increasing one-dimensional wavelength grid."""

    wave = np.asarray(wavelengths, dtype=float)
    if wave.ndim != 1:
        raise ValueError("wavelengths must be one-dimensional.")
    if wave.size < 2:
        raise ValueError("wavelengths must contain at least two points.")
    if not np.all(np.isfinite(wave)):
        raise ValueError("wavelengths must contain only finite values.")
    if np.any(np.diff(wave) <= 0.0):
        raise ValueError("wavelengths must be strictly increasing.")
    return wave


def validate_spectrum_grid(wavelengths, spectrum):
    """Validate that a spectrum uses `wavelengths` as its last dimension."""

    spectrum = np.asarray(spectrum)
    if spectrum.shape == ():
        raise ValueError("spectrum must have a wavelength dimension.")
    if spectrum.shape[-1] != len(wavelengths):
        raise ValueError("spectrum last dimension must match wavelengths.")


def _validate_positive_fwhm(fwhm):
    fwhm = float(fwhm)
    if not np.isfinite(fwhm) or fwhm <= 0.0:
        raise ValueError("Gaussian PSF FWHM must be positive.")
    return fwhm
