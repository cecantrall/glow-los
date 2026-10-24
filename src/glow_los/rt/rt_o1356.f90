module rt_o1356

use rt_utils

implicit none

contains

subroutine rt(ns,ds,bkss,vers,tablepath,toa,source,cont)

  implicit none
  
  ! Define I/O
  integer, intent(in) :: ns
  real, intent(in) :: ds(ns)
  real, intent(in) :: bkss(ns,3)
  real, intent(in) :: vers(ns)
  character (len=200), intent(in) :: tablepath
  real, intent(out) :: toa
  real, intent(out) :: source(ns), cont(ns)

  ! Grid variables
  real :: scdo2(ns) ! O2 slant column density above each path point
  real :: scdo(ns) ! O slant column density above each path point
  real :: init(ns)
  real :: xoxt(ns)
  real :: sssf(ns), dssf(ns), tssf(ns)
  real :: tsf1356(ns), tsf1359(ns) 
  real :: tn(ns), o(ns), o2(ns)

! O2 and O FUV cross sections
  integer, parameter :: bin_size = 1200
  real :: o2xsect_wave(bin_size), o2xsect(bin_size)
  real :: xox1356, xox1359
  real :: xo2_1356, xo2_1359
  real :: tau_O2_1356(ns), tau_O2_1359(ns), tau(ns)
  real :: tau_O_1356(ns), tau_O_1359(ns)
  real :: scb1356, scb1359
  real :: weight_1356(ns), weight_1359(ns)

! Define parameters
  real, parameter :: a1356 = 4.540E3 ! Transition coefficient for 135.6
  real, parameter :: a1359 = 1.340E3 ! Transition coefficient for 135.9
  real, parameter :: r1356 = a1356/(a1356+a1359)
  real, parameter :: r1359 = a1359/(a1356+a1359)
  real, parameter :: xox0 = 2.5e-18
  real, parameter :: xo2 = 7.0E-18 ! O2 absorption cross section, cm -2
  integer :: k, kk

  !! Define background variables
  o = bkss(:,1)
  o2 = bkss(:,2)
  tn = bkss(:,3)

!! Compute O2 optical depth
  call o2_scd(scdo2,o2,ds,ns)
  ! Read O2 absorption cross section data
  call read_o2fuvxsect(o2xsect_wave,o2xsect,tablepath)
  ! Interpolate O2 data on wavelength grid
  call interpolate_value(o2xsect_wave,o2xsect,1356.0,xo2_1356)
  call interpolate_value(o2xsect_wave,o2xsect,1359.0,xo2_1359)
  ! Compute optical depth at specific wavelengths
  tau_O2_1356 = scdo2*xo2_1356
  tau_O2_1359 = scdo2*xo2_1359

  !! Compute O optical depth
  ! Calculate O slant column density above each path point
  call o_scd(scdo,o,ds,ns)
  ! Calculate O cross sections
  call oxsect(xox0,xox1356,xox1359)
  ! Compute optical depth
  tau_O_1356 = scdo*xox1356
  tau_O_1359 = scdo*xox1359

  ! Apply temperature dependent ground state population to cross section
  xoxt = xox0*(0.812+(tn-200)*(0.617-0.812)/(1000.0-200.0))*1.15
  
  ! Calculate dz and optical depth in each layer
  tau = (o2*xo2+o*xoxt)*ds

  ! single scatttering of 135.6
  sssf=0.0
  init = vers*r1356
  do k = 2, ns-1
    do kk = 1, k
      sssf(k) = sssf(k) + xoxt(k)*o(k)*init(kk)*ds(kk)*exp(-sum(tau(kk:k)))
    end do
    do kk = k+2, ns
      sssf(k) = sssf(k) + xoxt(k)*o(k)*init(kk)*ds(kk)*exp(-sum(tau(k:kk)))
    end do
  end do
  sssf = sssf*r1356

  ! double scattering of 135.6
  dssf = 0.0
  do k = 2, ns-1
    do kk = 1, k
      dssf(k) = dssf(k) + xoxt(k)*o(k)*sssf(kk)*ds(kk)*exp(-sum(tau(kk:k)))
    end do
    do kk = k+2, ns
      dssf(k) = dssf(k) + xoxt(k)*o(k)*sssf(kk)*ds(kk)*exp(-sum(tau(k:kk)))
    end do
  end do
  dssf = dssf*r1356

  ! triple scattering of 135.6
  tssf=0.0
  do k = 2, ns-1
    do kk = 1, k
      tssf(k) = tssf(k) + xoxt(k)*o(k)*dssf(kk)*ds(kk)*exp(-sum(tau(kk:k)))
    end do
    do kk = k+2, ns
      tssf(k) = tssf(k) + xoxt(k)*o(k)*dssf(kk)*ds(kk)*exp(-sum(tau(k:kk)))
    end do
  end do
  tssf = tssf*r1356

  ! Add scattered photons proporionally to 135.6 and 135.9 total source functions:
  tsf1356 = init + sssf + dssf + tssf
  tsf1359 = vers*r1359 + (sssf+dssf+tssf)*r1359/r1356

  !! Calculate emergent intensity for 135.6 and 135.9
  do k = 1, ns
    weight_1356(k) = tsf1356(k)*ds(k)*exp(-tau_O2_1356(k)-tau_O_1356(k))
    weight_1359(k) = tsf1359(k)*ds(k)*exp(-tau_O2_1359(k)-tau_O_1359(k)) !mR to R
  end do
  source = tsf1356+tsf1359

  scb1356 =  sum(weight_1356)
  scb1359 =  sum(weight_1359)

  ! Total doublet intensity
  toa = scb1356+scb1359

  ! Total contribution
  cont = weight_1356+weight_1359

end subroutine rt

end module rt_o1356
