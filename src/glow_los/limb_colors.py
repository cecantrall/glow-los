from __future__ import annotations

import numpy as np

from . import PACKAGE_ROOT
from .extensions import ni, vk
from .los import los as LOSModel
from .spectral import (
    combine_emissions,
    convolve_spectrum,
    convolved_emission_spectrum,
    integrate_bands,
    to_rayleighs,
)


FUV_COLOR_BANDS = {
    "1356": (1343.0, 1377.0),
    "lbhs": (1410.0, 1528.0),
    "lbhl": (1671.0, 1812.0),
}

FUV_COLOR_BANDS_POST_2004 = {
    **FUV_COLOR_BANDS,
    "lbhl": (1672.0, 1792.0),
}


def update_glowparams(
    run_object,
    f107,
    f107a,
    f107p,
    iscale,
    xuvfac,
    kchem,
    jlocal,
    itail,
    fmono,
    emono,
    ef,
    ec,
    lbh_scale_fac,
) -> None:
    run_object.d_glow.update(
        {
            "f107": f107,
            "f107a": f107a,
            "f107p": f107p,
            "iscale": iscale,
            "xuvfac": xuvfac,
            "kchem": kchem,
            "jlocal": jlocal,
            "itail": itail,
            "fmono": fmono,
            "emono": emono,
            "ef": ef,
            "ec": ec,
            "lbh_scale_fac": lbh_scale_fac,
        }
    )


def simulate_limb_colors(
    lat,
    lon,
    sza,
    tanalts,
    bks,
    basepath=None,
    f107=120,
    f107a=120,
    f107p=120,
    iscale=1,
    xuvfac=3,
    kchem=2,
    jlocal=0,
    itail=0,
    fmono=0,
    emono=0,
    ef=0,
    ec=0,
    lbh_scale_fac=1.351,
    psf_fwhm=None,
    include_o1356_rr=False,
    include_vk=False,
    include_ni=False,
    color_bands=None,
    return_components=False,
):
    if basepath is None:
        basepath = str(PACKAGE_ROOT)

    simulator = LOSModel(basepath)
    update_glowparams(
        simulator,
        f107,
        f107a,
        f107p,
        iscale,
        xuvfac,
        kchem,
        jlocal,
        itail,
        fmono,
        emono,
        ef,
        ec,
        lbh_scale_fac,
    )

    _, ceta, eta = simulator.glow1d_only(lat, lon, sza, bks)
    tanalts = np.asarray(tanalts, dtype=float)
    selected_bands = FUV_COLOR_BANDS if color_bands is None else dict(color_bands)
    colors = {name: [None] * len(tanalts) for name in selected_bands}
    colors["alt"] = [None] * len(tanalts)
    if include_o1356_rr:
        if "1356" not in selected_bands:
            raise ValueError("include_o1356_rr requires a band named '1356'.")
        colors["1356_rr"] = [None] * len(tanalts)
    component_names = []
    if include_vk:
        component_names.append("vk")
    if include_ni:
        component_names.append("ni1743")
    if return_components:
        colors["components"] = {
            component: {
                name: [None] * len(tanalts) for name in selected_bands
            }
            for component in component_names
        }

    redister_cache = simulator.redister_o1356_cache_from_eta(bks, eta)

    wave = (
        simulator.d_rt["lbh_swave"]
        + simulator.d_rt["lbh_apb"]
        * np.arange(simulator.d_rt["lbh_nwave"])
    )
    extension_vers = {}
    if include_vk:
        extension_vers["vk"] = vk.from_glow(ceta, wave)
    if include_ni:
        extension_vers["ni1743"] = ni.from_glow(ceta, wave)

    for i, tanalt in enumerate(tanalts):
        toa1356, _, _ = simulator.los_1d_o1356_from_redister_cache(bks, *redister_cache, tan=tanalt)
        toalbh, _, _, native_wave = simulator.los_1d_from_eta("n2lbh", lat, lon, sza, bks, eta, ceta, tan=tanalt)
        toa1493, _, _ = simulator.los_1d_from_eta("n1493", lat, lon, sza, bks, eta, ceta, tan=tanalt)

        if not np.allclose(native_wave, wave):
            raise RuntimeError("Native LBH and extension wavelength grids differ.")
        extension_toa = {}
        for component, volume_emission_rate in extension_vers.items():
            extension_toa[component], _, _ = simulator.los_1d_from_spectral_ver(
                bks, wave, volume_emission_rate, tanalt
            )
        total_spectral_toa = toalbh + sum(
            extension_toa.values(), start=np.zeros_like(toalbh)
        )

        line_radiances = {1356.0: toa1356, 1493.0: toa1493}
        if psf_fwhm is None:
            spec = combine_emissions(wave, total_spectral_toa, line_radiances)
        else:
            spec = convolved_emission_spectrum(
                wave, total_spectral_toa, line_radiances, psf_fwhm
            )
        color_values = integrate_bands(wave, spec, selected_bands)

        for name, value in color_values.items():
            colors[name][i] = value
        colors["alt"][i] = tanalt
        if return_components:
            for component, toa in extension_toa.items():
                component_spectrum = to_rayleighs(toa)
                if psf_fwhm is not None:
                    component_spectrum = convolve_spectrum(
                        wave, component_spectrum, psf_fwhm
                    )
                component_values = integrate_bands(
                    wave, component_spectrum, selected_bands
                )
                for name, value in component_values.items():
                    colors["components"][component][name][i] = value
        if include_o1356_rr:
            toa1356_rr, _, _ = simulator.los_1d_from_eta("o1356_rr", lat, lon, sza, bks, eta, ceta, tan=tanalt)
            colors["1356_rr"][i] = color_values["1356"] * (
                toa1356_rr / toa1356
            )

    return colors
