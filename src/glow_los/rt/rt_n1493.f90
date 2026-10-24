module rt_n1493

use rt_utils

implicit none

contains

subroutine rt(ns,ds,bkss,vers,tablepath,toa,source,cont)

  implicit none
  
  ! Define I/O
  integer, intent(in) :: ns
  real, intent(in) :: ds(ns)
  real, intent(in) :: bkss(ns)
  real, intent(in) :: vers(ns)
  character (len=200), intent(in) :: tablepath
  real, intent(out) :: toa
  real, intent(out) :: source(ns)
  real, intent(out) :: cont(ns)

  ! Grid variables
  real :: scdo2(ns) ! O2 slant column density above each path point
  real :: o2(ns)
  

! O2 and O FUV cross sections
  integer, parameter :: bin_size = 1200
  real :: o2xsect_wave(bin_size), o2xsect(bin_size)
  real :: xo2_1493
  real :: tau_O2_1493(ns)
  real :: scb1493, cell_tau, foreground_tau, average_transmission
  integer :: k

  !! Define background variables
  o2 = bkss(:)

!! Compute O2 optical depth
  call o2_scd(scdo2,o2,ds,ns)
  ! Read O2 absorption cross section data
  call read_o2fuvxsect(o2xsect_wave,o2xsect,tablepath)
  ! Interpolate O2 data on wavelength grid
  call interpolate_value(o2xsect_wave,o2xsect,1493.0,xo2_1493)
  ! Compute optical depth at specific wavelengths
  tau_O2_1493 = scdo2*xo2_1493

  !! Compute TOA radiance for each wave bin
  do k = 1, ns
      cell_tau = o2(k)*ds(k)*xo2_1493
      foreground_tau = tau_O2_1493(k)-cell_tau
      if (abs(cell_tau) > 1.0e-6) then
        average_transmission = (1.0-exp(-cell_tau))/cell_tau
      else
        average_transmission = 1.0-cell_tau/2.0
      end if
      cont(k) = vers(k)*ds(k)*average_transmission*exp(-foreground_tau)
  end do

  source = vers
  toa = sum(cont)

end subroutine rt

end module rt_n1493
