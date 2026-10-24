# Emissions derived from GLOW output

`glow_los.extensions` contains optional emission calculations that consume
native GLOW output without modifying GLOW itself. They return altitude-resolved
volume-emission spectra in `cm^-3 s^-1` per wavelength bin. The generic
`los.los_1d_from_spectral_ver` method then applies the package O2 absorption
table and the same exact within-cell limb integration used by the native FUV
transport routines.

## Vegard-Kaplan bands

`glow_los.extensions.vk.from_glow(ceta, wavelengths)` expands GLOW's native
VK(0,9) component into an N2 A-X band spectrum. Band centers and Einstein A
coefficients come from the Gilmore, Laher, and Espy (1992) compilation bundled
in `data/n2_vk_bands.csv`.

GLOW does not provide vibrationally resolved A-state populations. The default
uses a normalized geometric population over v'=0..10, constrained so that
VK(0,9) exactly reproduces GLOW's native effective rate. This approximation
does not model cascade production or altitude-dependent quenching by O, N2,
O2, and NO. A caller can pass a different normalized population vector through
`vibrational_populations` without changing the LOS code.

## Atomic nitrogen 174.3 nm

`glow_los.extensions.ni.from_glow(ceta, wavelengths)` retains GLOW's separate
photodissociation and electron-impact N I 149.3 nm production components and
scales each to 174.3 nm. The default 174.3/149.3 ratio is `0.43/1.20`, from the
200 eV emission cross sections reported by Mumma and Zipf (1971). Separate
ratios can be supplied for the two production mechanisms.

This is an optically thin source calculation with O2 absorption along the LOS.
It does not model atomic-nitrogen resonant self-absorption.

## Limb-color integration

The existing behavior remains the default. Enable either source explicitly:

```python
from glow_los.limb_colors import (
    FUV_COLOR_BANDS_POST_2004,
    simulate_limb_colors,
)

colors = simulate_limb_colors(
    lat,
    lon,
    sza,
    tangent_altitudes,
    background,
    include_vk=True,
    include_ni=True,
    color_bands=FUV_COLOR_BANDS_POST_2004,
    return_components=True,
)
```

`return_components=True` adds separately integrated `vk` and `ni1743`
diagnostics. These are useful for auditing the contribution to any supplied
spectral band definition; the extension modules are not tied to GUVI or LBHL.

## References

- Gilmore, F. R., Laher, R. R., and Espy, P. J. (1992), *J. Phys. Chem. Ref.
  Data*, 21, 1005-1107.
- Mumma, M. J., and Zipf, E. C. (1971), *J. Chem. Phys.*, 55, 1661.
- Bishop, J., and Feldman, P. D. (2003), *J. Geophys. Res.*, 108(A6), 1243.
