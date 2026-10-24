"""Optional emissions derived from native GLOW outputs.

These modules do not modify GLOW.  They translate resolved GLOW production
components into additional volume-emission spectra that can be passed through
the standard :mod:`glow_los.los` radiative-transfer path.
"""

from . import ni, vk

__all__ = ["ni", "vk"]
