# Python wrapper for redister
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path

import numpy as np

MAX_REDISTER_LAYERS = 99


def redister_o1356_1d(
    bkss,
    vers,
    redister_src,
    redister_out="redister.out",
    return_weights=False,
):
    n_layers = len(vers)
    if n_layers > MAX_REDISTER_LAYERS:
        raise ValueError(
            f"REDISTER supports at most {MAX_REDISTER_LAYERS} layers; got {n_layers}. "
            "Use a coarser or shorter background altitude grid."
        )
    if np.asarray(bkss).shape[0] != n_layers:
        raise ValueError("bkss and vers must have the same number of altitude levels.")

    redister_src = Path(redister_src)
    with _redister_runtime_dir() as runtime_dir:
        # setup input files
        emission = "o1356"
        setup(emission, bkss, vers, redister_src, runtime_dir)

        # run
        process = run(redister_src, runtime_dir)
        returncode = process.wait()
        if returncode != 0:
            raise RuntimeError(f"REDISTER failed with exit code {returncode}.")
    
        # read output
        fsource, kext, alt, spectral_weights = read(
            runtime_dir / redister_out, len(vers)
        )
    
        # get source by summing over wavelengths
    #fsource_sum = np.sum(fsource,axis=1) # sum over spectral regions
    
    if return_weights:
        return fsource, kext, alt, spectral_weights
    return fsource, kext, alt


@contextmanager
def _redister_runtime_dir():
    root = _default_redister_runtime_root()
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="redister-", dir=root) as runtime_dir:
        yield Path(runtime_dir)


def _default_redister_runtime_root():
    return Path(tempfile.gettempdir()) / "glow-los"


def setup(emission, bkss, vers, redister_src, runtime_dir):
    redister_src = Path(redister_src)
    runtime_dir = Path(runtime_dir)
    # setup redister.inp (different template for different emission)
    if emission == "o1356":
        with open(redister_src / "redister_o1356.inp", "r") as file_template:
            Lines = file_template.readlines()
        
        # write to input file with modifications
        with open(runtime_dir / "redister.inp", "w+") as file_write:
            count = 0
            for line in Lines:
                if count == 0:
                    content_temp = line.split(",")
                    content_temp[0] = str(len(vers))
                    content_temp[1] = str(len(vers))
                    content_temp[2] = str(len(vers))
                    content = ",".join(map(str, content_temp))
                    file_write.write(content)
                else:
                    content = line.strip()
                    file_write.write(content + "\n")
                count += 1
            
        # setup redister.pro (1356 volume emission rates (ph/cm3/s)—aka the initial source function)
        #print(vers)
        with open(runtime_dir / "redister.pro", "w+") as file:
            for i in range(len(vers)):
                content = str(vers[i])
                file.write(content + "\n")

        # setup redister.atm (neutral atmosphere (cm-3), temperature (K) and altitude grid (km))
        # [alt, temp, [O], [N2], [O2], [N], [H]]
        with open(runtime_dir / "redister.atm", "w+") as file:
            for i in range(len(vers)):
                content = " ".join(map(str, bkss[i, :]))
                file.write(content + "\n")
    

def run(redister_src, runtime_dir):
    bashCommand = str(Path(redister_src) / "a.out")
    process = subprocess.Popen([bashCommand], cwd=runtime_dir)
    return process
    

def read(outfile, nlev):
    outfile = Path(outfile)
    # read lines from output file
    with open(outfile, "r") as file_template:
        Lines = file_template.readlines()

    spectral_weights = _read_spectral_weights(Lines, outfile)

    ## READ KEXT ##
    # Find start of kext line
    count=0
    kextstart = None
    for line in Lines:
        count+=1
        if line.startswith('   (N)    Z(CM)                KEXT(CM-1) ((M =  1 TO M = 12) L =  1 TO L =  2)'):
            kextstart = count
            continue
    # Find end of kext line
    count=0
    kextend = None
    for line in Lines:
        count+=1
        if line.startswith('   (N)    Z(CM)                TAU ((M =  1 TO M = 12) L =  1 TO L =  2)'):
            kextend = count-4
            continue
    if kextstart is None or kextend is None:
        raise ValueError(f"Could not parse KEXT section in {outfile}.")
    # Read kext lines
    count=kextstart
    a=0
    kext = np.empty((nlev,24))
    alt = np.empty((nlev))
    while count < kextend:
        if len(Lines[count].split()) > 6:
            if a >= nlev:
                break
            lst = Lines[count].split()[2:]+Lines[count+1].split()[:]+Lines[count+2].split()[:]+Lines[count+3].split()[:]
            kext[a,:] = [float(i) for i in lst]
            alt[a] = float(Lines[count].split()[1])
            count+=4
            a+=1
        else:
            count+=1
            
    ## READ final source ##
    # Find start of fsource line
    count=0
    fsourcestart = None
    for line in Lines:
        count+=1
        if line.startswith('   (N)    Z(CM)                FINAL SOURCE (((X =  1 TO X = 12) L =  1 TO L =  2) U = 1 TO U =  2)'):
            fsourcestart = count
            continue
    # Find end of fsource line
    count=0
    fsourceend = None
    for line in Lines:
        count+=1
        if line.startswith(' TOTAL FLUX OUT OF ATMOSPHERE (UP & DOWN)'):
            fsourceend = count-3
            continue
    if fsourcestart is None or fsourceend is None:
        raise ValueError(f"Could not parse FINAL SOURCE section in {outfile}.")
    # Read kext lines
    count=fsourcestart
    a=0
    fsource = np.empty((nlev,24))
    while count < fsourceend:
        if len(Lines[count].split()) > 6:
            if a >= nlev:
                break
            lst = Lines[count].split()[2:]+Lines[count+1].split()[:]+Lines[count+2].split()[:]+Lines[count+3].split()[:]
            fsource[a,:] = [float(i) for i in lst]
            count+=4
            a+=1
        else:
            count+=1
   
    return (
        np.flip(fsource, axis=0),
        np.flip(kext, axis=0),
        np.flip(alt),
        spectral_weights,
    )


def _read_spectral_weights(lines, outfile):
    """Read the 24 wavelength quadrature weights for the first angle ordinate."""

    start = next(
        (
            i + 1
            for i, line in enumerate(lines)
            if "MU WT(J)" in line and "X WT(J)" in line
        ),
        None,
    )
    if start is None:
        raise ValueError(f"Could not parse spectral weights in {outfile}.")

    weights = []
    for line in lines[start:]:
        fields = line.split()
        if len(fields) != 7:
            if weights:
                break
            continue
        try:
            row = int(fields[0])
            weight = float(fields[4])
        except ValueError:
            if weights:
                break
            continue
        if row > 24:
            break
        weights.append(weight)

    if len(weights) != 24:
        raise ValueError(
            f"Expected 24 REDISTER spectral weights in {outfile}; got {len(weights)}."
        )
    return np.asarray(weights)

## Functions that need work

# Not sure if this works as of 11/05/2024
#def redister_o1356_2d(ds,bkss,vers,alt_path,redister_src,redister_out):
#    # change to redister source path
#    cwd = os.getcwd()
#    os.chdir(redister_src)

#    # run redister looping over locs
#    fsource_path = np.empty((len(ds),24)); kext_path = np.empty((len(ds),24))
#    for i in range(len(ds)):
#        # setup input files
#        emission = 'o1356'
#        setup(emission,bkss[i,:,:],vers[:,i],redister_src)

#        # run
#        run(redister_src)
#        os.wait()
    
#        # read output
#        fsource, kext, alt = read(redister_src+redister_out,len(vers[:,i]))
        
#        # interpolate to path      
#        # log interpolation for fsource
#        f = inter.interp1d(alt*1e-5,np.log(fsource),kind='linear',axis=0)
#        fsource_path[i,:] = np.exp(f(alt_path[i]))
#        # log interpolation for kext
#        f = inter.interp1d(alt*1e-5,np.log(kext),kind='linear',axis=0)
#        kext_path[i,:] = np.exp(f(alt_path[i]))
        
#    fsource_sum = np.sum(fsource_path,axis=1) # sum over spectral regions

#    # integrate along los
#    toa = integrate(ds,fsource_path,kext_path)
    
#    os.chdir(cwd)
    
#    return toa, fsource_sum
       

# Below are routines taken from redister to use for the line of sight integration 
## These routines are currently not being used and need to eventually be tested to improve line-os-sight-integration

#def integs(x, y, z, n):
#    #EVALUATE RESULT = INTEGRAL FROM X(1) TO X(N) OF Y * Z, ASSUMING THAT Y
#    #VARIES EXPONENTIALLY WITH X, WHILE Z IS LINEAR IN X.
#    result = 0.0
#    for i in range(1, n):
#        dx = x[i] - x[i - 1]
#        dy = y[i] - y[i - 1]
#        za = z[i - 1]
#        zb = z[i]
#        ea = math.exp(-y[i - 1])
#        eb = math.exp(-y[i])
#        if dy >= 0.1 and za > 1.0e-20 and zb > 1.0e-20:
#            result += dx * (za * ea - zb * eb - (za - zb) * (ea - eb) / dy) / dy
#        else:
#            result += dx * (za * ea + zb * eb) / 2.0
#    return result

#def integ(x, y):
#    #EVALUATE RESULT = INTEGRAL FROM X(1) TO X(N) OF Y,  USING A CUBIC
#    #SPLINE FIT TO Y, PLUS THE SPLINE QUADRATURE FORMULA OF FORSYTHE,
#    #MALCOLM, AND MOLER, P. 90.
#    n = len(x)
#    b, c, d = spline(n, x, y)
#    result = 0.0
#    for i in range(1, n):
#        dx = x[i] - x[i - 1]
#        result += dx * ((y[i] + y[i - 1]) / 2.0 - dx**2 * (c[i] + c[i - 1]) / 12.0)
#    return result

#def spline(N, X, Y):
#    # CUBIC SPLINE ROUTINE FROM FORSYTHE, MALCOLM, AND MOLER, P. 77.
#    B = np.empty(N)
#    C = np.empty(N)
#    D = np.empty(N)
#    NM1 = N - 1
#    if N >= 2:
#        if N >= 3:
#            D[0] = X[1] - X[0]
#            C[1] = (Y[1] - Y[0]) / D[0]
#            for i in range(1, NM1):
#                D[i] = X[i + 1] - X[i]
#                B[i] = 2.0 * (D[i - 1] + D[i])
#                C[i + 1] = (Y[i + 1] - Y[i]) / D[i]
#                C[i] = C[i + 1] - C[i]
#            B[0] = -D[0]
#            B[N - 1] = -D[NM1 - 1]
#            C[0] = 0.0
#            C[N - 1] = 0.0
#            if N > 3:
#                C[0] = C[2] / (X[3] - X[1]) - C[1] / (X[2] - X[0])
#                C[N - 1] = C[N - 2] / (X[N - 1] - X[N - 3]) - C[N - 3] / (X[N - 2] - X[N - 4])
#                C[0] = C[0] * D[0] ** 2 / (X[3] - X[0])
#                C[N - 1] = -C[N - 1] * D[NM1 - 1] ** 2 / (X[N - 1] - X[N - 4])
#            for i in range(1, N):
#                T = D[i - 1] / B[i - 1]
#                B[i] = B[i] - T * D[i - 1]
#                C[i] = C[i] - T * C[i - 1]
#            C[N - 1] = C[N - 1] / B[N - 1]
#            for i in range(NM1 - 1, -1, -1):
#                C[i] = (C[i] - D[i] * C[i + 1]) / B[i]
#            B[N - 1] = (Y[N - 1] - Y[NM1 - 1]) / D[NM1 - 1] + D[NM1 - 1] * (C[NM1 - 1] + 2.0 * C[N - 1])
#            for i in range(NM1):
#                B[i] = (Y[i + 1] - Y[i]) / D[i] - D[i] * (C[i + 1] + 2.0 * C[i])
#                D[i] = (C[i + 1] - C[i]) / D[i]
#                C[i] = 3.0 * C[i]
#            C[N - 1] = 3.0 * C[N - 1]
#            D[N - 1] = D[NM1 - 1]
#        else:
#            B[0] = (Y[1] - Y[0]) / (X[1] - X[0])
#            C[0] = 0.0
#            D[0] = 0.0
#            B[1] = B[0]
#            C[1] = 0
#            D[1] = 0
            
#    return B, C, D
            
