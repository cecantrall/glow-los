import numpy as np
from datetime import datetime

def datetime_to_idate(dt: datetime) -> int:
    """
    Convert a Python datetime to the Fortran-style IDATE (YYYYDDD).

    Parameters
    ----------
    dt : datetime
        The datetime object to convert.

    Returns
    -------
    idate : int
        The date as an integer in YYYYDDD format.
    """
    year = dt.year
    day_of_year = dt.timetuple().tm_yday
    idate = year * 1000 + day_of_year
    return idate

def suncor(idate: int, ut: float):
    """
    Compute solar declination (sdec), right ascension (srasn), and Greenwich sidereal time (gst)
    in radians for a given date and UT seconds.

    Parameters
    ----------
    idate : int
        Date in format YYDDD or YYYYDDD (e.g., 2024315 for 2024, day 315).
    ut : float
        Universal Time in seconds.

    Returns
    -------
    sdec : float
        Solar declination in radians.
    srasn : float
        Solar right ascension in radians.
    gst : float
        Greenwich sidereal time in radians.
    """
    pi = np.pi
    fday = ut / 86400.0
    iyr = idate // 1000
    iday = idate - iyr * 1000

    # Handle Y2K fix: supports 1950–2050
    if iyr >= 1900:
        iyr -= 1900
    if iyr < 50:
        iyr += 100

    dj = 365 * iyr + (iyr - 1) // 4 + iday + fday - 0.5
    t = dj / 36525.0

    vl = np.mod(279.696678 + 0.9856473354 * dj, 360.0)
    gst = np.mod(279.696678 + 0.9856473354 * dj + 360.0 * fday + 180.0, 360.0) * pi / 180.0
    g = np.mod(358.475845 + 0.985600267 * dj, 360.0) * pi / 180.0

    slong = vl + (1.91946 - 0.004789 * t) * np.sin(g) + 0.020094 * np.sin(2.0 * g)
    obliq = (23.45229 - 0.0130125 * t) * pi / 180.0
    slp = (slong - 0.005686) * pi / 180.0

    sind = np.sin(obliq) * np.sin(slp)
    cosd = np.sqrt(1.0 - sind**2)

    sdec = np.arctan(sind / cosd)
    srasn = pi - np.arctan2(1.0 / np.tan(obliq) * sind / cosd, -np.cos(slp) / cosd)

    return sdec, srasn, gst


def solzen(idate: int, ut: float, glat: float, glong: float) -> float:
    """
    Compute Solar Zenith Angle (SZA) in degrees.

    Parameters
    ----------
    idate : int
        Date in format YYDDD or YYYYDDD (e.g., 2024315 for 2024, day 315).
    ut : float
        Universal time in seconds.
    glat : float
        Geographic latitude in degrees.
    glong : float
        Geographic longitude in degrees.

    Returns
    -------
    sza : float
        Solar zenith angle in degrees.
    """
    pi = np.pi
    rlat = np.radians(glat)
    rlong = np.radians(glong)

    sdec, srasn, gst = suncor(idate, ut)

    rh = srasn - (gst + rlong)
    cossza = np.sin(sdec) * np.sin(rlat) + np.cos(sdec) * np.cos(rlat) * np.cos(rh)
    cossza = np.clip(cossza, -1.0, 1.0)
    sza = np.degrees(np.arccos(cossza))

    return sza