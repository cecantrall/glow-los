import numpy as np

from glow_los.extensions import ni, vk
from glow_los.los import los


def test_ni_1743_scales_native_glow_components_separately():
    ceta = np.zeros((3, 14, 2))
    ceta[0, 13] = [10.0, 20.0]
    ceta[1, 13] = [3.0, 4.0]

    components = ni.components_from_glow(
        ceta,
        photodissociation_ratio=0.5,
        electron_impact_ratio=0.25,
    )

    np.testing.assert_allclose(components["photodissociation"], [5.0, 10.0])
    np.testing.assert_allclose(components["electron_impact"], [0.75, 1.0])


def test_ni_1743_spectrum_conserves_the_line_volume_emission_rate():
    ceta = np.zeros((3, 14, 2))
    ceta[0, 13] = [10.0, 20.0]
    ceta[1, 13] = [2.0, 4.0]
    wave = np.arange(1700.0, 1800.0, 0.4)

    spectrum = ni.from_glow(ceta, wave)

    expected = (ceta[0, 13] + ceta[1, 13]) * ni.NI_1743_TO_1493_RATIO
    np.testing.assert_allclose(spectrum.sum(axis=1), expected)


def test_vk_spectrum_reproduces_the_native_vk_09_anchor():
    ceta = np.zeros((3, 14, 2))
    ceta[1, 0] = [7.0, 11.0]
    wave = np.array([3352.5, 3353.0, 3353.5])

    spectrum = vk.from_glow(ceta, wave)

    np.testing.assert_allclose(spectrum[:, 1], ceta[1, 0])
    np.testing.assert_allclose(spectrum.sum(axis=1), ceta[1, 0])


def test_vk_default_populations_are_normalized_and_match_glow_anchor_rate():
    populations = vk.default_vibrational_populations()
    anchor = vk.bands()[
        (vk.bands()["v_upper"] == 0) & (vk.bands()["v_lower"] == 9)
    ][0]

    np.testing.assert_allclose(populations.sum(), 1.0)
    np.testing.assert_allclose(
        populations[0] * anchor["einstein_a_s"],
        vk.GLOW_VK_09_EFFECTIVE_RATE_S,
    )


def test_absorbing_spectrum_reduces_to_path_sum_without_absorption():
    ds = np.array([2.0, 3.0])
    source = np.array([[1.0, 10.0], [2.0, 20.0]])

    toa, contribution = los.integrate_absorbing_spectrum(
        ds,
        source,
        absorber_density=np.zeros(2),
        cross_sections=np.array([0.5, 1.0]),
    )

    expected = source * ds[:, None]
    np.testing.assert_allclose(contribution, expected)
    np.testing.assert_allclose(toa, expected.sum(axis=0))


def test_absorbing_spectrum_uses_exact_attenuation_inside_a_cell():
    toa, contribution = los.integrate_absorbing_spectrum(
        ds=np.array([2.0]),
        source=np.array([[4.0]]),
        absorber_density=np.array([1.0]),
        cross_sections=np.array([0.5]),
    )

    expected = 4.0 / 0.5 * (1.0 - np.exp(-1.0))
    np.testing.assert_allclose(contribution, [[expected]])
    np.testing.assert_allclose(toa, [expected])


def test_packaged_o2_absorption_table_interpolates_to_spectral_grid():
    model = los()

    cross_sections = model.o2_absorption_cross_sections(
        np.array([1320.0, 1493.0, 1743.0, 1839.6])
    )

    assert cross_sections.shape == (4,)
    assert np.all(np.isfinite(cross_sections))
    assert np.all(cross_sections >= 0.0)
