# line-of-sight calculations for 1-d and 2-d case
from glow_los.background import bkmodel
from glow_los.glow import glowpy
from glow_los.rt import rt_fuv
from glow_los.rt import rt_vis
from glow_los.rt import redister
import numpy as np
from netCDF4 import Dataset
from scipy import interpolate as inter
from . import PACKAGE_ROOT

class los:
    def __init__(self,basepath=None):
        if basepath is None:
            basepath = str(PACKAGE_ROOT)
        self.basepath = basepath
        
        # Define structure to store glow input parameters
        self.d_glow = {'f107': 120,\
            'f107a': 120,\
            'f107p': 120,\
            'iscale': 1,\
            'xuvfac': 3,\
            'kchem': 2,\
            'jlocal':0,\
            'itail':0,\
            'fmono':0,\
            'emono':0,\
            'ef':0,\
            'ec':0,\
            'data_dir':basepath+'/glow/data/',\
            'lbh_scale_fac':1.351}
        # Scaling factor comes from Correira et al. 2021 (used for GOLD algorithm)

        self.d_rt = {'tables':basepath+'/rt/tables/',\
            'redister_src':basepath+'/rt/redister_src/',\
            'redister_outfile':'redister.out',\
            'lbh_nwave':1300,\
            'lbh_apb':0.4,\
            'lbh_swave':1320.0}
        self._o2_absorption_table = None

    def glow1d_only(self,lat,lon,sza,bks):
        nlev = np.size(bks,axis=0)
        auroral = np.array([self.d_glow['itail'],self.d_glow['fmono'],self.d_glow['emono'],self.d_glow['ef'],self.d_glow['ec']],order='F')
        switches = np.array([self.d_glow['iscale'],self.d_glow['xuvfac'],self.d_glow['kchem'],self.d_glow['jlocal']],order='F')
        solar_geomags = np.array([self.d_glow['f107'],self.d_glow['f107a'],self.d_glow['f107p']],order='F')
        xden,ceta,eta = glowpy.glowdriver.glow1d(lat,lon,sza,nlev,bks, \
                                                switches,solar_geomags,auroral,self.d_glow['data_dir'])
        return xden,ceta,eta

    @staticmethod
    def redister_o1356_states(bks):
        # REDISTER expects [z,T,O,N2,O2,N,H]. The 14-column GLOW
        # background schema stores O+ in bks[:,11], not neutral H, so leave H
        # at zero until neutral H is carried as a separate column.
        neutral_h = np.zeros_like(bks[:,0])
        return np.vstack([bks[:,0],bks[:,4],bks[:,1],bks[:,3],bks[:,2],bks[:,10],neutral_h]).T

    def redister_o1356_cache_from_eta(self, bks, eta):
        """
        Run REDISTER once for the full O1356 vertical atmosphere.

        The returned source and extinction arrays are sorted on an ascending
        altitude grid in km for reuse by tangent-altitude LOS integrations.
        """

        bkss = self.redister_o1356_states(bks)
        source, kext, alt, spectral_weights = redister.redister_o1356_1d(
            bkss,
            eta[12, :],
            self.d_rt["redister_src"],
            self.d_rt["redister_outfile"],
            return_weights=True,
        )
        alt_km = self._redister_alt_to_km(alt)
        order = np.argsort(alt_km)
        return source[order, :], kext[order, :], alt_km[order], spectral_weights

    def los_1d_o1356_from_redister_cache(
        self,
        bks,
        redister_source,
        redister_kext,
        redister_alt_km,
        spectral_weights,
        tan,
    ):
        """Return O1356 limb radiance from a cached REDISTER source solution."""

        alts, source, kext = self._redister_cache_limb_los(
            redister_alt_km,
            redister_source,
            redister_kext,
            tan,
        )
        half_ds = self.limb_ds(alts, tan) / 2.0
        ds, source_path, kext_path = self._symmetric_limb_path(
            half_ds, source, kext
        )

        nlev = np.size(bks, axis=0)
        cont = np.full(nlev, np.nan)
        ss = nlev - len(half_ds)
        toa, path_cont = self.integrate_redister(
            np.asarray(ds) * 1e5,
            source_path,
            kext_path,
            spectral_weights,
        )
        cont_path = self._fold_symmetric_contribution(path_cont)
        if ss >= 0:
            cont[ss:] = cont_path
        else:
            cont[:] = cont_path[-nlev:]
        return toa, source, cont

    def los_1d_from_spectral_ver(self, bks, wavelengths, volume_emission_rate, tan):
        """Integrate an arbitrary emission spectrum along an O2-absorbing limb ray.

        ``volume_emission_rate`` must have shape ``(altitude, wavelength)`` and
        units ``cm^-3 s^-1`` per wavelength bin. The returned top-of-atmosphere
        spectrum is in the native GLOW/RT radiance units used elsewhere in this
        package; divide by ``1e6`` for Rayleighs.
        """

        bks = np.asarray(bks, dtype=float)
        wave = np.asarray(wavelengths, dtype=float)
        vers = np.asarray(volume_emission_rate, dtype=float)
        if bks.ndim != 2 or bks.shape[1] < 3:
            raise ValueError("bks must contain altitude and O2 density columns.")
        if wave.ndim != 1 or wave.size == 0 or not np.all(np.isfinite(wave)):
            raise ValueError("wavelengths must be a nonempty finite one-dimensional grid.")
        if vers.shape != (bks.shape[0], wave.size):
            raise ValueError("volume_emission_rate must have shape (altitude, wavelength).")

        alts, o2, source = self._limb_spectral_profiles(
            bks[:, 0], bks[:, 2], vers, tan
        )
        half_ds = self.limb_ds(alts, tan) / 2.0
        ds, o2_path, source_path = self._symmetric_limb_path(
            half_ds, o2, source
        )
        cross_sections = self.o2_absorption_cross_sections(wave)
        toa, path_cont = self.integrate_absorbing_spectrum(
            np.asarray(ds) * 1e5,
            source_path,
            o2_path,
            cross_sections,
        )
        return toa, source, self._fold_symmetric_contribution(path_cont)

    def _limb_spectral_profiles(self, alts, o2, source, tan):
        alts = np.asarray(alts, dtype=float)
        o2 = np.asarray(o2, dtype=float)
        source = np.asarray(source, dtype=float)
        if alts.size == 0:
            return alts, o2, source

        eps = 1e-6
        tan = min(max(float(tan), alts[0]), alts[-1] - eps)
        index = int(np.searchsorted(alts, tan, side="left"))
        if index <= 0:
            return alts, o2, source
        if index >= alts.size:
            index = alts.size - 1
        if np.isclose(alts[index], tan):
            return alts[index:], o2[index:], source[index:, :]

        alts_out = np.concatenate(([tan], alts[index:]))
        o2_out = np.concatenate(([np.interp(tan, alts, o2)], o2[index:]))
        source_tan = self._interp_columns(tan, alts, source)
        source_out = np.vstack([source_tan, source[index:, :]])
        return alts_out, o2_out, source_out

    def o2_absorption_cross_sections(self, wavelengths):
        """Interpolate the package O2 FUV cross sections onto a wavelength grid."""

        wave = np.asarray(wavelengths, dtype=float)
        if self._o2_absorption_table is None:
            filename = self.d_rt["tables"] + "o2absxsectfuv.nc"
            with Dataset(filename) as dataset:
                table_wave = np.asarray(
                    dataset.variables["WAVELENGTH"][:], dtype=float
                ).squeeze()
                table_cross_section = np.asarray(
                    dataset.variables["O2XSECT"][:], dtype=float
                ).squeeze()
            if table_wave.ndim != 1 or table_cross_section.shape != table_wave.shape:
                raise ValueError("O2 absorption table must contain matching vectors.")
            order = np.argsort(table_wave)
            self._o2_absorption_table = (
                table_wave[order],
                table_cross_section[order],
            )
        table_wave, table_cross_section = self._o2_absorption_table
        if np.any(wave < table_wave[0]) or np.any(wave > table_wave[-1]):
            raise ValueError("wavelengths fall outside the O2 absorption table.")
        return np.interp(wave, table_wave, table_cross_section)

    @staticmethod
    def _redister_alt_to_km(alt):
        alt = np.asarray(alt, dtype=float)
        if alt.size and np.nanmax(np.abs(alt)) > 1.0e4:
            return alt * 1.0e-5
        return alt

    def _redister_cache_limb_los(self, alts, source, kext, tan):
        alts = np.asarray(alts, dtype=float)
        source = np.asarray(source, dtype=float)
        kext = np.asarray(kext, dtype=float)
        tan = float(tan)
        if alts.size == 0:
            return alts, source, kext

        eps = 1e-6
        tan = min(max(tan, alts[0]), alts[-1] - eps)
        i = int(np.searchsorted(alts, tan, side="left"))
        if i <= 0:
            return alts, source, kext
        if i >= len(alts):
            i = len(alts) - 1

        if np.isclose(alts[i], tan):
            return alts[i:], source[i:, :], kext[i:, :]

        alts_out = np.concatenate(([tan], alts[i:]))
        source_tan = self._interp_columns(tan, alts, source)
        kext_tan = self._interp_columns(tan, alts, kext)
        source_out = np.vstack([source_tan, source[i:, :]])
        kext_out = np.vstack([kext_tan, kext[i:, :]])
        return alts_out, source_out, kext_out

    @staticmethod
    def _interp_columns(x, xp, values):
        values = np.asarray(values, dtype=float)
        if values.ndim == 1:
            return np.interp(x, xp, values)
        return np.array([np.interp(x, xp, values[:, j]) for j in range(values.shape[1])])

    def los_1d_from_eta(self,feature,lat,lon,sza,bks,eta,ceta,tan=None,oza=None,file=None,min_lev=None,timeind=None):
        """
        Like los_1d, but uses precomputed GLOW outputs (eta/ceta).

        Inputs
        feature - emission feature (o1356,o1493,lbh,red,green)
        time - datetime
        lat - latitudes
        lon - longitudes
        tan - tangent altitude (if limb vewing desired)
        oza - observing zenith angle (if disk viewing desired)
        bks - background states
        eta, ceta - GLOW outputs from glow1d
        file - filename (if want to read from file, otherwise MSIS is used)
        min_lev - if reading from file, must provide lowest level to read from
        timeind - if reading from file, must provide time index if more than one time step is in file (must correspond to same time in time input)
        
        Outputs
        toa
        source
        cont 
        alts
        """
        if (tan != None) and (oza != None):
            print('Only tan or oza can be specified')
            return
            
        elif (tan == None) and (oza == None):
            print('tan or oza needs to be specified')
            return

        nlev = np.size(bks,axis=0)

        # LOS integration
        if feature == 'o1356':
            bkss = self.redister_o1356_states(bks)
            
            if oza != None:
                # get alts, states, vers
                alts,states,vers = bks[:,0],bkss,eta[12,:]
                # compute ds
                ds = self.disk_ds(alts,oza)
            elif tan != None:
                # get alts, states, vers
                alts,states,vers = self.limb_los(bks[:,0],bkss,eta[12,:],tan)
                # compute ds
                ds = self.limb_ds(alts,tan)
            
            # los integration
            source, kext, alt, spectral_weights = redister.redister_o1356_1d(
                states,
                vers,
                self.d_rt['redister_src'],
                self.d_rt['redister_outfile'],
                return_weights=True,
            )

            if tan is not None:
                half_ds = np.asarray(ds) / 2.0
                ds, source_path, kext_path = self._symmetric_limb_path(
                    half_ds, source, kext
                )
            else:
                source_path, kext_path = source, kext
            
            # integrate along los
            cont = np.full(nlev,np.nan)
            ss = nlev-len(source)
            toa, cont_path = self.integrate_redister(
                np.array(ds) * 1e5,
                source_path,
                kext_path,
                spectral_weights,
            )
            if tan is not None:
                cont_path = self._fold_symmetric_contribution(cont_path)
            cont[ss:] = cont_path
            
        elif feature == 'o1356_rr':
            bkss = self.redister_o1356_states(bks)
            
            if oza != None:
                # get alts, states, vers
                alts,states,vers = bks[:,0],bkss,ceta[2,12,:]
                # compute ds
                ds = self.disk_ds(alts,oza)
            elif tan != None:
                # get alts, states, vers
                alts,states,vers = self.limb_los(bks[:,0],bkss,ceta[2,12,:],tan)
                # compute ds
                ds = self.limb_ds(alts,tan)
            
            # los integration
            source, kext, alt, spectral_weights = redister.redister_o1356_1d(
                states,
                vers,
                self.d_rt['redister_src'],
                self.d_rt['redister_outfile'],
                return_weights=True,
            )

            if tan is not None:
                half_ds = np.asarray(ds) / 2.0
                ds, source_path, kext_path = self._symmetric_limb_path(
                    half_ds, source, kext
                )
            else:
                source_path, kext_path = source, kext
            
            # integrate along los
            cont = np.full(nlev,np.nan)
            ss = nlev-len(source)
            toa, cont_path = self.integrate_redister(
                np.array(ds) * 1e5,
                source_path,
                kext_path,
                spectral_weights,
            )
            if tan is not None:
                cont_path = self._fold_symmetric_contribution(cont_path)
            cont[ss:] = cont_path
            
        elif feature == 'o1356_ss':
            # grab required states ([O,O2,T])
            bkss = np.vstack([bks[:,1], bks[:,2], bks[:,4]]).T
            
            if oza != None:
                # get alts, states, vers
                alts,states,vers = bks[:,0],bkss,eta[12,:]
                # compute ds
                ds = self.disk_ds(alts,oza)
            elif tan != None:
                # get alts, states, vers
                alts,states,vers = self.limb_los(bks[:,0],bkss,eta[12,:],tan)
                # compute ds
                ds = self.limb_ds(alts,tan)
            
            # los integration
            cont = np.full(nlev,np.nan)
            ss = nlev-len(ds)
            toa,source,cont[ss:] = rt_fuv.rt_o1356.rt(np.array(ds)*1e5,states,vers,self.d_rt['tables'])  
            
        elif feature == 'n2lbh':
            # grab required states ([O2,T])
            bkss = np.vstack([bks[:,2], bks[:,4]]).T
            
            if oza != None:
                # get alts, states, vers
                alts,states,vers = bks[:,0],bkss,eta[11,:]
                # compute ds
                ds = self.disk_ds(alts,oza)
            elif tan != None:
                # get alts, states, vers
                alts,states,vers = self.limb_los(bks[:,0],bkss,eta[11,:],tan)
                # compute ds
                ds = self.limb_ds(alts,tan)

            # define wavelegnths
            waves=[]
            for i in range(self.d_rt['lbh_nwave']):
                waves.append(self.d_rt['lbh_swave'] + self.d_rt['lbh_apb']*i)
            # scale volume emission rate
            vers = vers*self.d_glow['lbh_scale_fac']

            if tan is not None:
                half_ds = np.asarray(ds) / 2.0
                ds, states_path, vers_path = self._symmetric_limb_path(
                    half_ds, states, vers
                )
            else:
                states_path, vers_path = states, vers
            
            # los integration
            cont = np.full((nlev,len(waves)),np.nan)
            ss = nlev-len(vers)
            toa, source_path, cont_path = rt_fuv.rt_lbh.rt(
                np.array(ds) * 1e5,
                states_path,
                vers_path,
                waves,
                self.d_rt['tables'],
            )
            if tan is not None:
                source = source_path[len(vers):, :]
                cont_path = self._fold_symmetric_contribution(cont_path)
            else:
                source = source_path
            cont[ss:,:] = cont_path
        
        elif feature == 'n1493':
            # grab required states ([O2])
            bkss = np.array(bks[:,2])
            
            if oza != None:
                # get alts, states, vers
                alts,states,vers = bks[:,0],bkss,eta[13,:]
                # compute ds
                ds = self.disk_ds(alts,oza)
            elif tan != None:
                # get alts, states, vers
                alts,states,vers = self.limb_los(bks[:,0],bkss,eta[13,:],tan)
                # compute ds
                ds = self.limb_ds(alts,tan)

            if tan is not None:
                half_ds = np.asarray(ds) / 2.0
                ds, states_path, vers_path = self._symmetric_limb_path(
                    half_ds, states, vers
                )
            else:
                states_path, vers_path = states, vers

            # los integration
            cont = np.full(nlev,np.nan)
            ss = nlev-len(vers)
            toa, source_path, cont_path = rt_fuv.rt_n1493.rt(
                np.array(ds) * 1e5,
                states_path,
                vers_path,
                self.d_rt['tables'],
            )
            if tan is not None:
                source = source_path[len(vers):]
                cont_path = self._fold_symmetric_contribution(cont_path)
            else:
                source = source_path
            cont[ss:] = cont_path
        
        elif feature == 'o5577':
            ### Add code to do wind component projection ###
            
            # grab required states ([T,Vproj])
            bkss = np.vstack([bks[:,4],bks[:,13]]).T
            
            if oza != None:
                # get alts, states, vers
                alts,states,vers = bks[:,0],bkss,eta[3,:]
                # compute ds
                ds = self.disk_ds(alts,oza)
            elif tan != None:
                # get alts, states, vers
                alts,states,vers = self.limb_los(bks[:,0],bkss,eta[3,:],tan)
                # compute ds
                ds = self.limb_ds(alts,tan)
                
            # los integration
            toa,source,waves = rt_vis.o5577(np.array(ds)*1e5,states,vers)
        
        elif feature == 'o6300':   
            ### Add code to do wind component projection ###
            
            # grab required states ([T,Vproj])
            bkss = np.vstack([bks[:,4],bks[:,13]]).T
            
            if oza != None:
                # get alts, states, vers
                alts,states,vers = bks[:,0],bkss,eta[4,:]
                # compute ds
                ds = self.disk_ds(alts,oza)
            elif tan != None:
                # get alts, states, vers
                alts,states,vers = self.limb_los(bks[:,0],bkss,eta[4,:],tan)
                # compute ds
                ds = self.limb_ds(alts,tan)
                
            # los integration
            toa,source,waves = rt_vis.o6300(np.array(ds)*1e5,states,vers)
        
        else:
            print('Do not recognize feature')
        
        if (feature == 'o5577') or (feature == 'o6300') or (feature == 'n2lbh'):
            return toa, source, cont, waves
        else: 
            return toa, source, cont

    def los_1d(self,feature,lat,lon,sza,bks,tan=None,oza=None,file=None,min_lev=None,timeind=None):
        """
        Inputs
        feature - emission feature (o1356,o1493,lbh,red,green)
        time - datetime
        lat - latitudes
        lon - longitudes
        tan - tangent altitude (if limb vewing desired)
        oza - observing zenith angle (if disk viewing desired)
        bks - background states
        file - filename (if want to read from file, otherwise MSIS is used)
        min_lev - if reading from file, must provide lowest level to read from
        timeind - if reading from file, must provide time index if more than one time step is in file (must correspond to same time in time input)
        
        Outputs
        toa
        source
        cont 
        alts
        """ 
        # Check that either tan or oza is specified
        if (tan != None) and (oza != None):
            print('Only tan or oza can be specified')
            return
            
        elif (tan == None) and (oza == None):
            print('tan or oza needs to be specified')
            return
        
        xden,ceta,eta = self.glow1d_only(lat,lon,sza,bks)
        return self.los_1d_from_eta(feature,lat,lon,sza,bks,eta,ceta,tan=tan,oza=oza,file=file,min_lev=min_lev,timeind=timeind)


    ## Supporting functions
    def disk_ds(self,alts,oza):
        ds = np.empty(len(alts))
        for i in range(len(alts)-1):
            ds[i] = abs(alts[i+1] - alts[i])/np.cos(np.deg2rad(oza))
        ds[-1] = ds[-2]
        
        return ds      

    def limb_ds(self,alts,tan):
        Re = 6371.0 # radius of earth (km)
        ds = np.empty(len(alts))

        for i in range(len(alts)-1):
            r_t = Re + tan
            r1 = Re + alts[i]
            r2 = Re + alts[i+1]
            ds[i] = 2 * (np.sqrt(r2**2 - r_t**2) - np.sqrt(r1**2 - r_t**2))

        ds[-1] = ds[-2]
        return ds        

    def limb_los(self,alts,bkss,vers,tan):
        # Continuous LOS segment by interpolating at the tangent point
        alts = np.asarray(alts, dtype=float)
        tan = float(tan)
        if alts.size == 0:
            return alts, bkss, vers
        eps = 1e-6
        tan = min(max(tan, alts[0]), alts[-1] - eps)

        i = int(np.searchsorted(alts, tan, side="left"))
        if i <= 0:
            return alts, bkss, vers
        if i >= len(alts):
            i = len(alts) - 1

        if np.isclose(alts[i], tan):
            alts_out = alts[i:]
            if len(bkss.shape) > 1:
                bk_out = bkss[i:, :]
            else:
                bk_out = bkss[i:]
            vers_out = vers[i:]
            return alts_out, bk_out, vers_out

        alts_out = np.concatenate(([tan], alts[i:]))
        if len(bkss.shape) > 1:
            bk_tan = np.array([np.interp(tan, alts, bkss[:, j]) for j in range(bkss.shape[1])])
            bk_out = np.vstack([bk_tan, bkss[i:, :]])
        else:
            bk_tan = np.interp(tan, alts, bkss)
            bk_out = np.concatenate(([bk_tan], bkss[i:]))
        v_tan = np.interp(tan, alts, vers)
        vers_out = np.concatenate(([v_tan], vers[i:]))
        return alts_out, bk_out, vers_out

    @staticmethod
    def _symmetric_limb_path(half_ds, *profiles):
        """Order a symmetric limb ray from the far side to the observer."""

        half_ds = np.asarray(half_ds)
        path_ds = np.concatenate((half_ds[::-1], half_ds))
        path_profiles = []
        for profile in profiles:
            profile = np.asarray(profile)
            path_profiles.append(np.concatenate((profile[::-1], profile), axis=0))
        return (path_ds, *path_profiles)

    @staticmethod
    def _fold_symmetric_contribution(path_cont):
        """Combine far- and near-side path contributions by altitude level."""

        path_cont = np.asarray(path_cont)
        nlev = path_cont.shape[0] // 2
        return path_cont[:nlev][::-1] + path_cont[nlev:]
    
    def integrate_redister(self, ds, fsource, kext, spectral_weights=None):
        """Integrate an ordered path, including exact attenuation within each cell."""

        ds = np.asarray(ds, dtype=float)
        fsource = np.asarray(fsource, dtype=float)
        kext = np.asarray(kext, dtype=float)
        if fsource.shape != kext.shape or fsource.shape[0] != ds.size:
            raise ValueError("ds, fsource, and kext must describe the same path.")

        if spectral_weights is None:
            spectral_weights = np.ones(fsource.shape[1])
        spectral_weights = np.asarray(spectral_weights, dtype=float)
        if spectral_weights.shape != (fsource.shape[1],):
            raise ValueError("spectral_weights must contain one value per source column.")

        fcont = np.empty_like(fsource)
        foreground_transmission = np.ones(fsource.shape[1])
        for i in range(ds.size - 1, -1, -1):
            cell_tau = ds[i] * kext[i]
            average_transmission = np.ones_like(cell_tau)
            finite_tau = np.abs(cell_tau) > 1e-10
            average_transmission[finite_tau] = (
                -np.expm1(-cell_tau[finite_tau]) / cell_tau[finite_tau]
            )
            fcont[i] = (
                fsource[i]
                * ds[i]
                * average_transmission
                * foreground_transmission
            )
            foreground_transmission *= np.exp(-cell_tau)

        weighted_cont = fcont * spectral_weights
        return np.sum(weighted_cont), np.sum(weighted_cont, axis=1)

    @staticmethod
    def integrate_absorbing_spectrum(
        ds, source, absorber_density, cross_sections
    ):
        """Integrate a spectral source through one absorbing constituent."""

        ds = np.asarray(ds, dtype=float)
        source = np.asarray(source, dtype=float)
        absorber_density = np.asarray(absorber_density, dtype=float)
        cross_sections = np.asarray(cross_sections, dtype=float)
        if source.ndim != 2 or source.shape[0] != ds.size:
            raise ValueError("source must have shape (path, wavelength).")
        if absorber_density.shape != ds.shape:
            raise ValueError("absorber_density must contain one value per path cell.")
        if cross_sections.shape != (source.shape[1],):
            raise ValueError("cross_sections must contain one value per wavelength.")

        contribution = np.empty_like(source)
        foreground_transmission = np.ones(source.shape[1])
        for index in range(ds.size - 1, -1, -1):
            cell_tau = ds[index] * absorber_density[index] * cross_sections
            average_transmission = np.ones_like(cell_tau)
            finite_tau = np.abs(cell_tau) > 1e-10
            average_transmission[finite_tau] = (
                -np.expm1(-cell_tau[finite_tau]) / cell_tau[finite_tau]
            )
            contribution[index] = (
                source[index]
                * ds[index]
                * average_transmission
                * foreground_transmission
            )
            foreground_transmission *= np.exp(-cell_tau)
        return contribution.sum(axis=0), contribution
