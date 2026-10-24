#!/bin/bash

glow=.

if [[ -n ${NETCDFF_INCLUDE} ]]; then
  echo "NETCDFF_INCLUDE is set."
else
  echo "NETCDFF_INCLUDE is not set."
  echo "##########################################################################"
  echo
  echo "Please set NETCDFF_INCLUDE to the directory where netcdf.mod is located"
  echo
  exit
fi


if [[ -n ${NETCDFF_LIB} ]]; then
  echo "NETCDFF_LIB is set."
else
  echo "NETCDFF_LIB is not set."
  echo "##########################################################################"
  echo
  echo "Please set NETCDFF_LIB to the directory where libnetcdff.* is located"
  echo
  exit
fi
PYTHON_BIN=${PYTHON:-python}
export FC=${FC:-gfortran}


${PYTHON_BIN} -m numpy.f2py --build-dir $glow/build -c --opt='-O3 -ffree-line-length-none --std=legacy -O0 -fbacktrace' \
     -I${NETCDFF_INCLUDE} -L${NETCDFF_LIB} -lnetcdff  -m glowpy \
     $glow/cglow.f90 $glow/glowdriver.f90 $glow/glow.f90 $glow/bands.f90 \
     $glow/egrid.f90 $glow/ephoto.f90 $glow/etrans.f90 $glow/exsect.f $glow/fieldm.f $glow/gchem.f90 \
     $glow/maxt.f90 $glow/qback.f90 $glow/rcolum.f90  \
     $glow/ssflux.f90
