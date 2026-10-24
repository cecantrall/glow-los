def waccmx_1d(fname,lat,lon,min_lev,timeind=None):
    waccmx = nc.Dataset(fname)

    # get lat, lon lev
    lats = waccmx['lat'][:].data
    lons = waccmx['lon'][:].data
    P = np.squeeze(waccmx['lev'][:min_lev]).data

    # find lat and lon index
    latind = abs(lats-lat).argmin()
    lonind = abs(lons-lon).argmin() # lon should be between 0-360
    
    bks = np.zeros([min_lev,numvars],order='F')
    if timeind == None:
        try:
            Z3 = np.squeeze(waccmx['Z3'][:min_lev,latind,lonind]).data
            T = np.squeeze(waccmx['T'][:min_lev,latind,lonind]).data
            O2 = np.squeeze(waccmx['O2'][:min_lev,latind,lonind]).data
            O = np.squeeze(waccmx['O'][:min_lev,latind,lonind]).data
            N2 = 1-O-O2
            
            z = np.flip(geopotential_to_geometric(Z3))
            o = np.flip(molfraction_to_numberdens(O,P,T))
            n2 = np.flip(molfraction_to_numberdens(N2,P,T))
            o2 = np.flip(molfraction_to_numberdens(O2,P,T))
            
            bks[:,0] = z
            bks[:,1] = o # O
            bks[:,2] = o2 # O2
            bks[:,3] = n2 # N2
            bks[:,4] = np.flip(T) # T
        except:
            print('waccmx does not contain needed vars')
        # Less neccessary vars
        try:
            NE = np.squeeze(waccmx['EDens'][:min_lev,latind,lonind]).data
            bks[:,5] = np.flip(NE)
        except:
            print('NE does not exist - setting to 0')
        try:
            TE = np.squeeze(waccmx['TElec'][:min_lev,latind,lonind]).data
            bks[:,6] = np.flip(TE)
        except:
            print('TE does not exist - setting to 0')
        try:
            TI = np.squeeze(waccmx['TIon'][:min_lev,latind,lonind]).data
            bks[:,7] = np.flip(TI)
        except:
            print('TI does not exist - setting to 0')
        try:
            NO = np.squeeze(waccmx['NO'][:min_lev,latind,lonind]).data
            no = np.flip(molfraction_to_numberdens(NO,P,T))
            bks[:,8] = no
        except:
            print('NO does not exist - setting to 0')
        try:
            N = np.squeeze(waccmx['N'][:min_lev,latind,lonind]).data
            n4s = np.flip(molfraction_to_numberdens(N,P,T))
            bks[:,10] = n4s
        except:
            print('N does not exist - setting to 0')
        try:
            H = np.squeeze(waccmx['H'][:min_lev,latind,lonind]).data
            h = np.flip(molfraction_to_numberdens(H,P,T))
            bks[:,11] = h
        except:
            print('H does not exist - setting to 0')
        try:
            UN = np.squeeze(waccmx['U'][:min_lev,latind,lonind]).data
            bks[:,12] = np.flip(UN)
        except:
            print('UN does not exist - setting to 0')
        try:
            VN = np.squeeze(waccmx['V'][:min_lev,latind,lonind]).data
            bks[:,13] = np.flip(VN)
        except:
            print('VN does not exist - setting to 0')
            
    else:
        try:
            Z3 = np.squeeze(waccmx['Z3'][timeind,:min_lev,latind,lonind]).data
            T = np.squeeze(waccmx['T'][timeind,:min_lev,latind,lonind]).data
            O2 = np.squeeze(waccmx['O2'][timeind,:min_lev,latind,lonind]).data
            O = np.squeeze(waccmx['O'][timeind,:min_lev,latind,lonind]).data
            N2 = 1-O-O2
            
            z = np.flip(geopotential_to_geometric(Z3))
            o = np.flip(molfraction_to_numberdens(O,P,T))
            n2 = np.flip(molfraction_to_numberdens(N2,P,T))
            o2 = np.flip(molfraction_to_numberdens(O2,P,T))
            
            bks[:,0] = z
            bks[:,1] = o # O
            bks[:,2] = o2 # O2
            bks[:,3] = n2 # N2
            bks[:,4] = np.flip(T) # T
        except:
            print('waccmx does not contain needed vars')
        # Less neccessary vars
        try:
            NE = np.squeeze(waccmx['EDens'][timeind,:min_lev,latind,lonind]).data
            bks[:,5] = np.flip(NE)
        except:
            print('NE does not exist - setting to 0')
        try:
            TE = np.squeeze(waccmx['TElec'][timeind,:min_lev,latind,lonind]).data
            bks[:,6] = np.flip(TE)
        except:
            print('TE does not exist - setting to 0')
        try:
            TI = np.squeeze(waccmx['TIon'][timeind,:min_lev,latind,lonind]).data
            bks[:,7] = np.flip(TI)
        except:
            print('TI does not exist - setting to 0')
        try:
            NO = np.squeeze(waccmx['NO'][timeind,:min_lev,latind,lonind]).data
            no = np.flip(molfraction_to_numberdens(NO,P,T))
            bks[:,8] = no
        except:
            print('NO does not exist - setting to 0')    
        try:
            N = np.squeeze(waccmx['N'][timeind,:min_lev,latind,lonind]).data
            n4s = np.flip(molfraction_to_numberdens(N,P,T))
            bks[:,10] = n4s
        except:
            print('N does not exist - setting to 0')
        try:
            H = np.squeeze(waccmx['H'][timeind,:min_lev,latind,lonind]).data
            h = np.flip(molfraction_to_numberdens(H,P,T))
            bks[:,11] = h
        except:
            print('H does not exist - setting to 0')
        try:
            UN = np.squeeze(waccmx['U'][timeind,:min_lev,latind,lonind]).data
            bks[:,12] = np.flip(UN)
        except:
            print('UN does not exist - setting to 0')
        try:
            VN = np.squeeze(waccmx['V'][timeind,:min_lev,latind,lonind]).data
            bks[:,13] = np.flip(VN)
        except:
            print('VN does not exist - setting to 0')
        
                 
        return bks

def waccmx_2d(fname,lats,lons,min_lev,timeind=None):
    waccmx = nc.Dataset(fname)

    # get lat, lon lev
    latss = waccmx['lat'][:].data
    lonss = waccmx['lon'][:].data
    P = np.squeeze(waccmx['lev'][:min_lev]).data

    bks = np.zeros([len(lons),min_lev,numvars],order='F')
    for i in range(len(lons)):
        # find lat and lon index
        latind = abs(latss-lats[i]).argmin()
        lonind = abs(lonss-lons[i]).argmin() # lon should be between 0-360
    
        if timeind == None:
            try:
                Z3 = np.squeeze(waccmx['Z3'][:min_lev,latind,lonind]).data
                T = np.squeeze(waccmx['T'][:min_lev,latind,lonind]).data
                O2 = np.squeeze(waccmx['O2'][:min_lev,latind,lonind]).data
                O = np.squeeze(waccmx['O'][:min_lev,latind,lonind]).data
                N2 = 1-O-O2
            
                z = np.flip(geopotential_to_geometric(Z3))
                o = np.flip(molfraction_to_numberdens(O,P,T))
                n2 = np.flip(molfraction_to_numberdens(N2,P,T))
                o2 = np.flip(molfraction_to_numberdens(O2,P,T))
            
                bks[i,:,0] = z
                bks[i,:,1] = o # O
                bks[i,:,2] = o2 # O2
                bks[i,:,3] = n2 # N2
                bks[i,:,4] = np.flip(T) # T
            except:
                print('waccmx does not contain needed vars')
            # Less neccessary vars
            try:
                NE = np.squeeze(waccmx['EDens'][:min_lev,latind,lonind]).data
                bks[i,:,5] = np.flip(NE)
            except:
                print('NE does not exist - setting to 0')
            try:
                TE = np.squeeze(waccmx['TElec'][:min_lev,latind,lonind]).data
                bks[i,:,6] = np.flip(TE)
            except:
                print('TE does not exist - setting to 0')
            try:
                TI = np.squeeze(waccmx['TIon'][:min_lev,latind,lonind]).data
                bks[i,:,7] = np.flip(TI)
            except:
                print('TI does not exist - setting to 0')
            try:
                NO = np.squeeze(waccmx['NO'][:min_lev,latind,lonind]).data
                no = np.flip(molfraction_to_numberdens(NO,P,T))
                bks[i,:,8] = no
            except:
                print('NO does not exist - setting to 0')
            try:
                N = np.squeeze(waccmx['N'][:min_lev,latind,lonind]).data
                n4s = np.flip(molfraction_to_numberdens(N,P,T))
                bks[i,:,10] = n4s
            except:
                print('N does not exist - setting to 0')
            try:
                H = np.squeeze(waccmx['H'][:min_lev,latind,lonind]).data
                h = np.flip(molfraction_to_numberdens(H,P,T))
                bks[i,:,11] = h
            except:
                print('H does not exist - setting to 0')
            try:
                UN = np.squeeze(waccmx['U'][:min_lev,latind,lonind]).data
                bks[i,:,12] = np.flip(UN)
            except:
                print('UN does not exist - setting to 0')
            try:
                VN = np.squeeze(waccmx['V'][:min_lev,latind,lonind]).data
                bks[i,:,13] = np.flip(VN)
            except:
                print('VN does not exist - setting to 0')
            
        else:
            try:
                Z3 = np.squeeze(waccmx['Z3'][timeind,:min_lev,latind,lonind]).data
                T = np.squeeze(waccmx['T'][timeind,:min_lev,latind,lonind]).data
                O2 = np.squeeze(waccmx['O2'][timeind,:min_lev,latind,lonind]).data
                O = np.squeeze(waccmx['O'][timeind,:min_lev,latind,lonind]).data
                N2 = 1-O-O2
            
                z = np.flip(geopotential_to_geometric(Z3))
                o = np.flip(molfraction_to_numberdens(O,P,T))
                n2 = np.flip(molfraction_to_numberdens(N2,P,T))
                o2 = np.flip(molfraction_to_numberdens(O2,P,T))
            
                bks[i,:,0] = z
                bks[i,:,1] = o # O
                bks[i,:,2] = o2 # O2
                bks[i,:,3] = n2 # N2
                bks[i,:,4] = np.flip(T) # T
            except:
                print('waccmx does not contain needed vars')
            # Less neccessary vars
            try:
                NE = np.squeeze(waccmx['EDens'][timeind,:min_lev,latind,lonind]).data
                bks[i,:,5] = np.flip(NE)
            except:
                print('NE does not exist - setting to 0')
            try:
                TE = np.squeeze(waccmx['TElec'][timeind,:min_lev,latind,lonind]).data
                bks[i,:,6] = np.flip(TE)
            except:
                print('TE does not exist - setting to 0')
            try:
                TI = np.squeeze(waccmx['TIon'][timeind,:min_lev,latind,lonind]).data
                bks[i,:,7] = np.flip(TI)
            except:
                print('TI does not exist - setting to 0')
            try:
                NO = np.squeeze(waccmx['NO'][timeind,:min_lev,latind,lonind]).data
                no = np.flip(molfraction_to_numberdens(NO,P,T))
                bks[i,:,8] = no
            except:
                print('NO does not exist - setting to 0')    
            try:
                N = np.squeeze(waccmx['N'][timeind,:min_lev,latind,lonind]).data
                n4s = np.flip(molfraction_to_numberdens(N,P,T))
                bks[i,:,10] = n4s
            except:
                print('N does not exist - setting to 0')
            try:
                H = np.squeeze(waccmx['H'][timeind,:min_lev,latind,lonind]).data
                h = np.flip(molfraction_to_numberdens(H,P,T))
                bks[i,:,11] = h
            except:
                print('H does not exist - setting to 0')
            try:
                UN = np.squeeze(waccmx['U'][timeind,:min_lev,latind,lonind]).data
                bks[i,:,12] = np.flip(UN)
            except:
                print('UN does not exist - setting to 0')
            try:
                VN = np.squeeze(waccmx['V'][timeind,:min_lev,latind,lonind]).data
                bks[i,:,13] = np.flip(VN)
            except:
                print('VN does not exist - setting to 0')
            
                 
    return bks

def waccmx_3d(fname,min_lev,timeind=None):
    waccmx = nc.Dataset(fname)

    lat = waccmx['lat'][:].data
    lon = waccmx['lon'][:].data
    
    # get lev
    P = np.squeeze(waccmx['lev'][:min_lev]).data
    Pin = np.empty((min_lev,len(lat),len(lon)))
    for i in range(len(lat)):
        for j in range(len(lon)):
            Pin[:,i,j]=P

    bks = np.zeros([len(lat),len(lon),min_lev,numvars],order='F')
    if timeind == None:
        Z3 = np.squeeze(waccmx['Z3'][0,:min_lev,:,:]).data
        z = np.flip(geopotential_to_geometric(Z3),axis=0)
        
        try:
            Z3 = np.squeeze(waccmx['Z3'][0,:min_lev,:,:]).data
            T = np.squeeze(waccmx['T'][0,:min_lev,:,:]).data
            O2 = np.squeeze(waccmx['O2'][0,:min_lev,:,:]).data
            O = np.squeeze(waccmx['O'][0,:min_lev,:,:]).data
            N2 = 1-O-O2
        
            z = np.flip(geopotential_to_geometric(Z3),axis=0)
            o = np.flip(molfraction_to_numberdens(O,Pin,T),axis=0)
            n2 = np.flip(molfraction_to_numberdens(N2,Pin,T),axis=0)
            o2 = np.flip(molfraction_to_numberdens(O2,Pin,T),axis=0)
            
            bks[:,:,:,0] = np.einsum('ijk->jki',z)
            bks[:,:,:,1] = np.einsum('ijk->jki',o) # O
            bks[:,:,:,2] = np.einsum('ijk->jki',o2) # O2
            bks[:,:,:,3] = np.einsum('ijk->jki',n2) # N2
            bks[:,:,:,4] = np.einsum('ijk->jki',np.flip(T,axis=0)) # T
            
            
        except:
            print('waccmx does not contain needed vars')
        # Less neccessary vars
        try:
            NE = np.squeeze(waccmx['EDens'][0,:min_lev,:,:]).data
            bks[:,:,:,5] = np.einsum('ijk->jki',np.flip(NE,axis=0))
        except:
            print('NE does not exist - setting to 0')
        try:
            TE = np.squeeze(waccmx['TElec'][0,:min_lev,:,:]).data
            bks[:,:,:,6] = np.einsum('ijk->jki',np.flip(TE,axis=0))
        except:
            print('TE does not exist - setting to 0')
        try:
            TI = np.squeeze(waccmx['TIon'][0,:min_lev,:,:]).data
            bks[:,:,:,7] = np.einsum('ijk->jki',np.flip(TI,axis=0))
        except:
            print('TI does not exist - setting to 0')
        try:
            NO = np.squeeze(waccmx['NO'][0,:min_lev,:,:]).data
            no = np.flip(molfraction_to_numberdens(NO,Pin,T),axis=0)
            bks[:,:,:,8] = np.einsum('ijk->jki',no)
        except:
            print('NO does not exist - setting to 0')
        try:
            N = np.squeeze(waccmx['N'][0,:min_lev,:,:]).data
            n4s = np.flip(molfraction_to_numberdens(N,Pin,T),axis=0)
            bks[:,:,:,10] = np.einsum('ijk->jki',n4s)
        except:
            print('N does not exist - setting to 0')
        try:
            H = np.squeeze(waccmx['H'][0,:min_lev,:,:]).data
            h = np.flip(molfraction_to_numberdens(H,Pin,T),axis=0)
            bks[:,:,:,11] = np.einsum('ijk->jki',h)
        except:
            print('H does not exist - setting to 0')
        try:
            UN = np.squeeze(waccmx['U'][0,:min_lev,:,:]).data
            bks[:,:,:,12] = np.flip(UN)
        except:
            print('UN does not exist - setting to 0')
        try:
            VN = np.squeeze(waccmx['V'][0,:min_lev,:,:]).data
            bks[:,:,:,13] = np.flip(VN)
        except:
            print('VN does not exist - setting to 0')
          
    else:
        try:
            Z3 = np.squeeze(waccmx['Z3'][timeind,:min_lev,:,:]).data
            T = np.squeeze(waccmx['T'][timeind,:min_lev,:,:]).data
            O2 = np.squeeze(waccmx['O2'][timeind,:min_lev,:,:]).data
            O = np.squeeze(waccmx['O'][timeind,:min_lev,:,:]).data
            N2 = 1-O-O2
            
            z = np.flip(geopotential_to_geometric(Z3),axis=0)
            o = np.flip(molfraction_to_numberdens(O,Pin,T),axis=0)
            n2 = np.flip(molfraction_to_numberdens(N2,Pin,T),axis=0)
            o2 = np.flip(molfraction_to_numberdens(O2,Pin,T),axis=0)
            
            bks[:,:,:,0] = np.einsum('ijk->jki',z)
            bks[:,:,:,1] = np.einsum('ijk->jki',o) # O
            bks[:,:,:,2] = np.einsum('ijk->jki',o2) # O2
            bks[:,:,:,3] = np.einsum('ijk->jki',n2) # N2
            bks[:,:,:,4] = np.einsum('ijk->jki',np.flip(T,axis=0)) # T
        except:
            print('waccmx does not contain needed vars')
        # Less neccessary vars
        try:
            NE = np.squeeze(waccmx['EDens'][timeind,:min_lev,:,:]).data
            bks[:,:,:,5] = np.einsum('ijk->jki',np.flip(NE,axis=0))
        except:
            print('NE does not exist - setting to 0')
        try:
            TE = np.squeeze(waccmx['TElec'][timeind,:min_lev,:,:]).data
            bks[:,:,:,6] = np.einsum('ijk->jki',np.flip(TE,axis=0))
        except:
            print('TE does not exist - setting to 0')
        try:
            TI = np.squeeze(waccmx['TIon'][timeind,:min_lev,:,:]).data
            bks[:,:,:,7] = np.einsum('ijk->jki',np.flip(TI,axis=0))
        except:
            print('TI does not exist - setting to 0')
        try:
            NO = np.squeeze(waccmx['NO'][timeind,:min_lev,:,:]).data
            no = np.flip(molfraction_to_numberdens(NO,Pin,T),axis=0)
            bks[:,:,:,8] = np.einsum('ijk->jki',no)
        except:
            print('NO does not exist - setting to 0')   
        try:
            N = np.squeeze(waccmx['N'][timeind,:min_lev,:,:]).data
            n4s = np.flip(molfraction_to_numberdens(N,Pin,T),axis=0)
            bks[:,:,:,10] = np.einsum('ijk->jki',n4s)
        except:
            print('N does not exist - setting to 0')
        try:
            H = np.squeeze(waccmx['H'][timeind,:min_lev,:,:]).data
            h = np.flip(molfraction_to_numberdens(H,Pin,T),axis=0)
            bks[:,:,:,11] = np.einsum('ijk->jki',h)
        except:
            print('H does not exist - setting to 0')
        try:
            UN = np.squeeze(waccmx['U'][timeind,:min_lev,:,:]).data
            bks[:,:,:,12] = np.flip(UN)
        except:
            print('UN does not exist - setting to 0')
        try:
            VN = np.squeeze(waccmx['V'][timeind,:min_lev,:,:]).data
            bks[:,:,:,13] = np.flip(VN)
        except:
            print('VN does not exist - setting to 0')
                   
    return bks

# functions used for WACCMX processing
def geopotential_to_geometric(z3):
    er = 6356.75 # km
    
    z3 = z3*1e-3 # m to km
    z = (er*z3)/(er-z3) 
    
    return z

def molfraction_to_numberdens(molfrac,p,t):
    k = 1.38e-23
    n = 1e-6*molfrac*p*1e2/(k*t)
    
    return n