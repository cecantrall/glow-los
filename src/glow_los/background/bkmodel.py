from pymsis import msis
from pymsis import utils
from iricore import iri
import numpy as np
import netCDF4 as nc

def bks_from_msis_iri(time,lat,lon,levs,f107=None,f107a=None,ap=None):
    levs = np.asarray(levs, dtype=float)
    bks = np.zeros([len(levs),14],order='F')
    
    o,o2,n2,t,no,h = msis1d(time,lat,lon,levs,f107=f107,f107a=f107a,ap=ap)
    edens,te,ti,oplus,o2plus,noplus = iri1d_on_grid(time,lat,lon,levs,f107=f107,f107a=f107a)
    
    bks[:,0] = levs
    bks[:,1] = o
    bks[:,2] = o2
    bks[:,3] = n2
    bks[:,4] = t
    bks[:,5] = edens
    bks[:,6] = te
    bks[:,7] = ti
    bks[:,8] = no
    bks[:,9] = 0.
    bks[:,10] = 0.
    bks[:,11] = oplus
    bks[:,12] = o2plus
    bks[:,13] = noplus

    return bks


def iri1d(time,lat,lon,alts,f107=None,f107a=None):
    """
    Get IRI profile
    
    Inputs:
    time - datetime object
    alts - altitudes in km
    lat - latitude
    lon - longitude
    
    Outputs:
    edens - electron density cm-3
    te - electron temperature [K]
    ti - ion temperature [K]
    oplus - O+ [cm-3]
    o2plus - O2+ [cm-3]
    noplus - NO+ [cm-3]
    """
    
    df = iri(time,alts,lat,lon,version=20,**_iri_solar_kwargs(f107=f107,f107a=f107a))
    edens = df.edens*1e-6
    ti = df.itemp
    te = df.etemp

    # Treat IRI ion outputs as percentages (0-100).
    o_frac = np.clip(np.asarray(df.o, dtype=float) / 100.0, 0.0, 1.0)
    o2_frac = np.clip(np.asarray(df.o2, dtype=float) / 100.0, 0.0, 1.0)
    no_frac = np.clip(np.asarray(df.no, dtype=float) / 100.0, 0.0, 1.0)
    oplus = edens * o_frac
    o2plus = edens * o2_frac
    noplus = edens * no_frac
    
    return edens, te, ti, oplus, o2plus, noplus


def iri1d_on_grid(time,lat,lon,levs,f107=None,f107a=None):
    """
    Return IRI outputs interpolated onto the requested altitude grid.

    iricore accepts altitude input as [start, stop, step], not arbitrary
    nonuniform grids. Run IRI on a uniform support grid and resample onto levs.
    """

    levs = np.asarray(levs, dtype=float)
    if levs.size == 0:
        raise ValueError("levs must contain at least one altitude.")
    if levs.size == 1:
        altrange = [float(levs[0]), float(levs[0]), 1.0]
    else:
        diffs = np.diff(levs)
        if np.any(diffs <= 0.0):
            raise ValueError("levs must be strictly increasing.")
        altrange = [float(levs[0]), float(levs[-1]), float(np.min(diffs))]

    df = iri(time,altrange,lat,lon,version=20,**_iri_solar_kwargs(f107=f107,f107a=f107a))
    iri_alt = np.asarray(df.height, dtype=float)
    edens = _interp_profile(iri_alt, np.asarray(df.edens, dtype=float) * 1e-6, levs, log_space=True)
    ti = _interp_profile(iri_alt, np.asarray(df.itemp, dtype=float), levs, log_space=False)
    te = _interp_profile(iri_alt, np.asarray(df.etemp, dtype=float), levs, log_space=False)

    o_frac = np.clip(_interp_profile(iri_alt, np.asarray(df.o, dtype=float) / 100.0, levs, log_space=False), 0.0, 1.0)
    o2_frac = np.clip(_interp_profile(iri_alt, np.asarray(df.o2, dtype=float) / 100.0, levs, log_space=False), 0.0, 1.0)
    no_frac = np.clip(_interp_profile(iri_alt, np.asarray(df.no, dtype=float) / 100.0, levs, log_space=False), 0.0, 1.0)

    oplus = edens * o_frac
    o2plus = edens * o2_frac
    noplus = edens * no_frac
    return edens, te, ti, oplus, o2plus, noplus


def _iri_solar_kwargs(f107=None, f107a=None):
    kwargs = {}
    if f107 is not None:
        kwargs["oarr40"] = float(f107)
    if f107a is not None:
        kwargs["oarr45"] = float(f107a)
    elif f107 is not None:
        kwargs["oarr45"] = float(f107)
    return kwargs


def _interp_profile(z_src, y_src, z_dst, log_space=False):
    z_src = np.asarray(z_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    z_dst = np.asarray(z_dst, dtype=float)
    finite = np.isfinite(z_src) & np.isfinite(y_src)
    z_src = z_src[finite]
    y_src = y_src[finite]
    if z_src.size == 0:
        return np.full_like(z_dst, np.nan, dtype=float)
    if z_src.size == 1:
        return np.full_like(z_dst, y_src[0], dtype=float)

    order = np.argsort(z_src)
    z_src = z_src[order]
    y_src = y_src[order]

    if log_space:
        floor = 1e-30
        return np.exp(np.interp(z_dst, z_src, np.log(np.clip(y_src, floor, None))))
    return np.interp(z_dst, z_src, y_src)
    

def msis1d(time,lat,lon,levs,f107=None,f107a=None,ap=None):
    str_date = np.datetime64(time)

    if f107 == None:
        f107,f107a,ap = utils.get_f107_ap(str_date)
        temp = msis.run(str_date,lon,lat,levs,f107,f107a,[ap])
    else:
        ap = [[ap] * 7]
        temp = msis.run(str_date,lon,lat,levs,f107s=f107,f107as=f107a,aps=ap)

    temp = np.squeeze(temp)
    o = temp[:,3]*1e-6 # O
    o2 = temp[:,2]*1e-6 # O2
    n2 = temp[:,1]*1e-6 # N2
    t = temp[:,10] # T
    no = temp[:,9]*1e-6 # NO
    #n4s = 0 #temp[:,7]*1e-6 # N4S (ground state N)
    h = temp[:,5]*1e-6 # H
    
    return o, o2, n2, t, no, h

def save_bks_to_txt(filename, bks):
    # Define headers for readability
    headers = [
        "Altitude", "O", "O2", "N2", "T", "eDens", "Te", "Ti", 
        "NO", "N2D", "N4S", "O+", "O2+", "NO+"
    ]
    
    # Save as a tab-delimited text file with headers
    np.savetxt(
        filename,
        bks,
        header="\t".join(headers),
        fmt="%.6e",  # scientific notation with 6 decimal places
        delimiter="\t",
        comments=""  # no comment character at start of header
    )
    print(f"✅ BKS array saved to {filename}")


def read_bks_txt(filepath, as_dataframe=False):
    """
    Reads a BKS text file created by save_bks_to_txt().

    Parameters
    ----------
    filepath : str
        Path to the text file.
    as_dataframe : bool, optional
        If True, returns a pandas DataFrame instead of (array, headers).

    Returns
    -------
    np.ndarray or pd.DataFrame, list[str]
        If as_dataframe=False → (bks_array, headers)
        If as_dataframe=True → DataFrame with columns labeled.
    """
    import pandas as pd
    
    # Read the header line (first line starting with column names)
    with open(filepath, 'r') as f:
        first_line = f.readline().strip()

    # Split headers (tab-delimited)
    headers = first_line.replace("#", "").split()
    headers = [h.strip() for h in headers if h.strip()]

    # Load data (skip header row)
    bks = np.array(np.loadtxt(filepath, skiprows=1))

    if as_dataframe:
        df = pd.DataFrame(bks, columns=headers)
        return df
    else:
        return bks, headers
    
# Example usage:
if __name__ == "__main__":
    import numpy as np
    from utils.config import read_config
    from glow_los.background import bkmodel
    import sys
    
    # read config file
    config_path = sys.argv[1]
    cfg = read_config(config_path)
        
    # Define altitudes
    levs = np.arange(80, 510, 10)  # e.g., altitude levels in km

    # Compute the array
    bks = bkmodel.bks_from_msis_iri(cfg['time'], cfg['lat'], cfg['lon'], levs)

    # Save to text file
    bkmodel.save_bks_to_txt(cfg['bks_file'], bks)
