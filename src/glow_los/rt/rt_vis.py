# line of sight integration for red and green line
import numpy as np

def o6300(ds,bkss,vers):
    tn = bkss[:,0]
    v = bkss[:,1]
    
    # Define wavelength
    lam0 = 630.0 #nm
    lammin = 629.999
    lammax = 630.001
    lamres = 0.000001
    lams = np.arange(lammin,lammax,lamres)
    
    # mass of O
    M = 15.999
    N0 = 6.02e23
    m = M/N0
    
    cont = np.empty((len(lams),len(ds)))
    for i in range(len(ds)):
        # compute doppler broadening
        lshape = doppler_broadening(vers[i],tn[i],lam0*1e-9,lams*1e-9,m)
        
        # compute doppler shift
        v=0
        lamshift = doppler_shift(v,lam0*1e-9)
        binshift=lamshift/lamres
        lshape = np.roll(lshape, int(binshift))
        
        # integrate ofver line path
        cont[:,i] = lshape*ds[i]
    
    source = vers       
    toa = np.sum(cont,axis=1)  
    
    return toa,source,lams

def o5577(ds,bkss,vers):
    tn = bkss[:,0]
    v = bkss[:,1]
    
    # Define wavelength
    lam0 = 557.7 #nm
    lammin = 557.6999
    lammax = 557.7001
    lamres = 0.000001
    lams = np.arange(lammin,lammax,lamres)
    
    # mass of O
    M = 15.999
    N0 = 6.02e23
    m = M/N0
    
    cont = np.empty((len(lams),len(ds)))
    for i in range(len(ds)):
        # compute doppler broadening
        lshape = doppler_broadening(vers[i],tn[i],lam0*1e-9,lams*1e-9,m)
        
        # compute doppler shift
        v=0
        lamshift = doppler_shift(v,lam0*1e-9)
        binshift=lamshift/lamres
        lshape = np.roll(lshape, int(binshift))
        
        # integrate ofver line path
        cont[:,i] = lshape*ds[i]
    source = vers
    toa = np.sum(cont,axis=1)  
    
    return toa,source,lams


def doppler_shift(v,lam0):
    c = 2.99792e8 #m/s
    f0 =c/lam0
    
    f = (1+v/c)*f0
    lam = c/f
    return lam

def doppler_broadening(ver,tn,lam0,lams,m):
    c = 2.99792e8 #m/s
    kb = 1.381e-23 #J/K
    f0 = c/lam0
    
    fgrid = c/lams
    alpha = f0*np.sqrt(2*kb*tn/(m*c**2))
    
    lshape = np.empty(len(lams))
    for i in range(len(lams)):
        lshape[i] = (ver/(alpha*np.sqrt(np.pi)))*np.exp(-(fgrid[i]-f0)**2/alpha**2)
        
        
    return lshape