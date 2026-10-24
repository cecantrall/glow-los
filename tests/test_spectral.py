import numpy as np

from glow_los.spectral import convolved_emission_spectrum, gaussian_line_spectrum


def test_gaussian_line_spectrum_conserves_line_radiance() -> None:
    wave = np.arange(1300.0, 1400.5, 0.5)
    radiance = 1.2e6

    spectrum = gaussian_line_spectrum(wave, {1356.0: radiance}, fwhm=3.0)

    assert spectrum.shape == wave.shape
    assert np.isclose(np.sum(spectrum), radiance)
    assert np.isclose(wave[np.argmax(spectrum)], 1356.0)


def test_convolved_emission_spectrum_combines_lbh_and_lines_in_rayleighs() -> None:
    wave = np.arange(1300.0, 1400.5, 0.5)
    lbh = np.full_like(wave, 1.0e6)

    spectrum = convolved_emission_spectrum(wave, lbh, {1356.0: 2.0e6}, psf_fwhm=3.0)

    assert spectrum.shape == wave.shape
    assert np.all(np.isfinite(spectrum))
    assert np.isclose(np.sum(spectrum), wave.size + 2.0)
