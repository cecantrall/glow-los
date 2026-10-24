module rt_lbh

use rt_utils

implicit none

contains

subroutine rt(ns,ds,bkss,vers,nwave,waves,tablepath,toa,source,cont)

  implicit none
  
  ! Define I/O
  integer, intent(in) :: ns
  integer, intent(in) :: nwave
  real, intent(in) :: ds(ns)
  real, intent(in) :: bkss(ns,2)
  real, intent(in) :: vers(ns)
  real, intent(in) :: waves(nwave)
  character (len=200), intent(in) :: tablepath
  real, intent(out) :: toa(nwave)
  real, intent(out) :: source(ns,nwave)
  real, intent(out) :: cont(ns,nwave)

  ! lbhspec variables
  real :: lbhtable(3000,10,2000)
  real :: wavetable(3000)
  real :: lbh(10,ns)
  integer :: k,w
  real :: vplbh(10)
  real :: spectrum(nwave)
  real :: spec(ns,nwave)
  real :: apb, cell_tau, foreground_tau, average_transmission
  real :: wavelbh(nwave)

  ! Grid variables
  real :: scdo2(ns) ! O2 slant column density above each path point
  real :: tau_O2_lbh(ns,nwave) ! optical depth
  real :: tn(ns), o2(ns)

  ! O2 FUV cross sections
  integer, parameter :: bin_size = 1200
  real :: o2xsect_wave(bin_size), o2xsect(bin_size)
  real :: xo2(nwave) ! interpolated O2 FUV cross sections

  !! Define background variables
  o2 = bkss(:,1)
  tn = bkss(:,2)
 
  !! Read LBH lookup table
  call read_lbhtable(lbhtable,wavetable,tablepath)

  !! Spectral synthesis
  call bands(ns,vers,lbh)

  !! Compute spectrum at each level
  apb = (waves(nwave)-waves(1))/nwave
  do k = 1, ns
    vplbh = lbh(:,k)
    call lbhspec(lbhtable,wavetable,waves(1),waves(nwave),nwave,apb,tn(k),vplbh,wavelbh,spectrum)
    spec(k,:) = spectrum
  end do

  !! Absorption
  call o2_scd(scdo2,o2,ds,ns)
  ! Read O2 absorption cross section data
  call read_o2fuvxsect(o2xsect_wave,o2xsect,tablepath)
  ! Interpolate O2 data on wavelength grid
  call interpolate_1D(o2xsect_wave,o2xsect,waves,xo2)
  ! Compute optical depth
  do k = 1, ns
    tau_O2_lbh(k,:) = scdo2(k)*xo2
  end do
  
  !! Compute TOA radiance for each wave bin
  do w = 1, nwave
    do k = 1, ns
      cell_tau = o2(k)*ds(k)*xo2(w)
      foreground_tau = tau_O2_lbh(k,w)-cell_tau
      if (abs(cell_tau) > 1.0e-6) then
        average_transmission = (1.0-exp(-cell_tau))/cell_tau
      else
        average_transmission = 1.0-cell_tau/2.0
      end if
      cont(k,w) = spec(k,w)*ds(k)*average_transmission*exp(-foreground_tau)
    end do
  end do

  source = spec
  toa = sum(cont,DIM=1)

end subroutine rt

! -------------------------------------------------------------------------------

subroutine read_lbhtable(lbhtable, wavelength, lbhtable_file)
  use netcdf
  
  ! This is the name of the data file to be read
  character (len = 200), intent(in) :: lbhtable_file
  
  ! Reading 1D data for two parameters
  real, intent(out) :: lbhtable(3000,10,2000), wavelength(3000)
  
  ! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid, istat

  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file.
  istat = nf90_open(trim(lbhtable_file)//'lbhtable.nc', NF90_NOWRITE, ncid)
  
  ! Get the varid of the wavelength variable, based on its name.
  istat = nf90_inq_varid(ncid, "LBHTABLE", varid)
  ! Read the data.
  istat = nf90_get_var(ncid, varid, lbhtable) 
  
  ! Get the varid of the wavelength variable, based on its name.
  istat = nf90_inq_varid(ncid, "WAVELENGTH", varid) 
  ! Read the data.
  istat = nf90_get_var(ncid, varid, wavelength) 
  
  ! Close the file, freeing all resources.
  istat = nf90_close(ncid)
    
end subroutine read_lbhtable

!---------------------------------------------------------------------------

subroutine lbhspec(lbhtable, wavetable, wavemin, wavemax, nwave, apb, t, vplbh, wave, spect)
  ! INPUTS:
  !       wavemin:  minimum wavelength, integer, A (must be between 1200 and 2399)
  !       wavemax:  maximum wavelength, integer, A (must be between 1200 and 2399)
  !       t:        temperature, K, real
  !       vplbh:    array of lbh v'=0-9 population rates, cm-3 s-1, real, 10 elements
  ! OUTPUTS:  
  !       wave:     wavelengths of spectrum, 0.4 A bins, each bin centered at midpoints, real
  !       spect:    spectrum in 1-A bins, cm-3 s-1 per bin, real
  !
  ! PROCEDURE:
  !       Reads LBH lookup table from file 'lbhtable.nc' finds spectrum at given temperature for
  !       each upper state v'=0-9 of LBH bands, sums output spectrum for all bands.
  !       (Note that v'=7-9 excitation rates are generally zero because they predissociate). 
  !       Each v' spectra at each t in the table is normalized so that the sum of the bins = 1.0
  !       User is expected to provide v' populations in appropriate ratios
  ! REFERENCES:
  !       LBH model written by Scott Budzien
  ! MODIFICATION HISTORY:
  !       Stan Solomon, 3/15, earlier version using output from Cleary model
  !       Stan Solomon, 1/16, no changes except to lbhtable.nc, now based on Budzien model
  !       Stan Solomon, 8/16, corrected comment on units of output spect
  !       Stan Solomon, 2/17, changed from 1 A to 0.4 A resolution
  implicit none
  
  ! Define inputs and outputs
  integer, intent(in) :: nwave
  real, intent(in) :: apb
  real, intent(in) :: lbhtable(:,:,:), wavetable(:), wavemin, wavemax, t
  real, intent(in) :: vplbh(:)
  real, intent(out) :: wave(:), spect(:)

  ! Define internal variables
  integer :: w1, w2, tf
  
  ! Specify wavelength range and temperature of interest (1 is added for consistency with IDL code. Not sure if lbhtable is designed with IDL in mind)
  w1 = floor((wavemin-1200.19)/apb) + 1
  w2 = floor((wavemax-1200.19)/apb) + 1
  tf = floor(t) + 1

  ! Create output arrays
  wave = wavetable(w1:w2)
  spect = matmul(lbhtable(w1:w2,:,tf),vplbh) 
  
end subroutine lbhspec

! -------------------------------------------------------------------------------------

subroutine bands(ns,ver,zlbh)
  ! Subroutine BANDS

  ! This software is part of the GLOW model.  Use is governed by the Open Source
  ! Academic Research License Agreement contained in the file glowlicense.txt.
  ! For more information see the file glow.txt.

  ! Stan Solomon, 3/2015

  ! Version 0.1 is really just a stub that divides total N2 singlet state
  ! band system excitation rate (a, a', w) into 7 upper states of the N2
  ! LBH band system (v' = 0-6).  States above v'=6 are presumed to dissociate.
  ! Frank Condon factors from Ajello and Shemansky, JGR, 90, 9845-9861, 1985

  ! Input:
  !   Use-associated variables from cglow:
  !     jmax   number of altitude levels
  !     nc     number of emission components (in this case, v' levels)
  !     aglw   excited state array (state, species, altitude)
  !     zeta   airglow volume emission rate array (emission, species, altitude)
  ! Output:
  !   Use-associated variables from cglow:
    !     zlbh   total excitation rate to each LBH v' level (cm-3 s-1)


  implicit none

  integer :: k,m
  integer, parameter :: nc=10
  real :: vp_pop(nc)
  data vp_pop /0.068,0.132,0.168,0.18,0.158,0.105,0.07,0.0,0.0,0.0/ ! Determined from GOLD (Ajello et al. 2020)

  integer, intent(in) :: ns
  real, intent(in) :: ver(ns)
  real, intent(out) :: zlbh(nc,ns)

  do k=1,ns
    do m=1,nc
      zlbh(m,k)=ver(k)*vp_pop(m)
    end do
  end do

  return

end subroutine bands

end module rt_lbh
