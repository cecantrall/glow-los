import numpy as np

from glow_los.los import los
from glow_los.rt.redister import _read_spectral_weights


def test_redister_spectral_weights_are_read_from_first_angle_ordinate():
    lines = [
        "   (J)    MU(J)       X(J)     MU WT(J)   X WT(J)      F(J)       B(J)\n"
    ]
    lines.extend(
        f" {row:4d}  2.113E-01  0.000E+00  5.000E-01  {row:.3E}  0.000E+00  7.554E-01\n"
        for row in range(1, 49)
    )

    weights = _read_spectral_weights(lines, "redister.out")

    np.testing.assert_allclose(weights, np.arange(1.0, 25.0))


def test_integrate_redister_applies_spectral_weights():
    model = los.__new__(los)
    ds = np.array([2.0, 3.0])
    source = np.array([[1.0, 10.0], [2.0, 20.0]])
    kext = np.zeros_like(source)
    weights = np.array([0.25, 0.75])

    toa, contribution = model.integrate_redister(ds, source, kext, weights)

    expected = np.sum(source * ds[:, None] * weights, axis=1)
    np.testing.assert_allclose(contribution, expected)
    np.testing.assert_allclose(toa, np.sum(expected))


def test_integrate_redister_uses_exact_cell_attenuation():
    model = los.__new__(los)
    ds = np.array([2.0])
    source = np.array([[4.0]])
    kext = np.array([[0.5]])

    toa, contribution = model.integrate_redister(ds, source, kext)

    expected = source[0, 0] / kext[0, 0] * (1.0 - np.exp(-1.0))
    np.testing.assert_allclose(contribution, [expected])
    np.testing.assert_allclose(toa, expected)


def test_symmetric_limb_path_can_be_folded_to_altitude_grid():
    half_ds = np.array([1.0, 2.0, 3.0])
    profile = np.array([10.0, 20.0, 30.0])

    path_ds, path_profile = los._symmetric_limb_path(half_ds, profile)

    np.testing.assert_allclose(path_ds, [3.0, 2.0, 1.0, 1.0, 2.0, 3.0])
    np.testing.assert_allclose(path_profile, [30.0, 20.0, 10.0, 10.0, 20.0, 30.0])
    np.testing.assert_allclose(
        los._fold_symmetric_contribution(np.arange(1.0, 7.0)),
        [7.0, 7.0, 7.0],
    )
