#!/bin/bash
rt=.


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

${PYTHON_BIN} -m numpy.f2py -c --opt='-O3 -ffree-line-length-none -std=legacy' \
   -I${NETCDFF_INCLUDE} -L${NETCDFF_LIB} -lnetcdff  \
   -m rt_fuv $rt/rt_utils.f90 $rt/rt_lbh.f90 $rt/rt_o1356.f90 $rt/rt_n1493.f90
 

 
