module glowdriver

! Inputs
! idate_in - date in yyddd format
! ut_in - ut in seconds
! glat_in - latitude
! glong_in - longitude (-180 to 180)
! nlev - number of altitude/pressure levels
! bks - background variables (11 variables)
!   - zg (must be provided)
!   - o (must be provided)
!   - o2 (must be provided)
!   - n2 (must be provided)
!   - tn (must be provided)
!   - ne
!   - te 
!   - ti
!   - no
!   - n2d (currently always set to zero)
!   - n4s 
!   - h
!   - un
!   - vn
! switches (4)
!   1. iscale_in - 0 = Hinteregger solar model, 1 = EUVAC solar model, 2 = user specified
!   2. xuvfac_in - solar flux multiplicative factor, EUVAC 18-50 A or Hinteregger 18-250 A
!   3. kchem_in - 3 = use specified Ne, 4 = calculate Ne < 200 km, specified Ne > 200 km
!   With the exception of KCHEM=0 (no calculations), none are foolproof.
!        One irritation that may be noted when KCHEM=4 (calculate Ne below 200 but use
!        provided Ne above 200) is a discontinuity in Ne(z) at 200 km.  GCHEM attempts
!        to deal with this by interpolating from the 200 km level to three grid
!        points above it.  Another caveat with KCHEM=4 or KCHEM=3 (Ne provided,
!        ions calculated) is that if there is an incompatibility between the specifed
!        Ne(z) and the calculated ionization rates, negative values for some ions may
!        result.  KCHEM=2 (electrons and major ions provided, minor ions calculated) and
!        KCHEM=1 (electrons and all ions except O+(2D,2P) provided) should be fairly
!        reliable.
!   4. jlocal_in - 0 = full electron transport, 1 = local only
! solar_geomags (4)
!   1. f107_in - f107
!   2. f107a_in - f107a
!   3. f107p_in - f107p_in 
! auroral (5)
!   1. itail_in - Low-energy tail (1=yes, 0=no) 
!   2. fmono_in - Monoenergetic energy flux, erg cm-2
!   3. emono_in - Monoenergetic characteristic energy, keV
!   4. ef_in - Auroal electron energy flux, erg cm-2
!   5. ec_in - Auroral electron characteristic energy, eV


use cglow

implicit none

contains

subroutine glow1d(glat_in,glong_in,sza_in,nlev,bks,switches,solar_geomags,auroral,data_path,xden,ceta,eta)

implicit none

! Inputs
!*integer, intent(in) :: idate_in
real, intent(in) :: glong_in,glat_in,sza_in
real, intent(in) :: solar_geomags(3)
integer, intent(in) :: switches(4)
real, intent(in) :: auroral(5)
integer, intent(in) :: nlev
character(len=1024), intent(in) :: data_path
!f2py depend(nlev) bks, xden, ceta, eta, zg
integer, parameter :: nbks = 14
!f2py depend(nbks) bks
real, intent(in) :: bks(nlev,nbks)

!!Outputs
real, intent(out) :: xden(nex,nlev)
real, intent(out) :: ceta(nc,nw,nlev)
real, intent(out) :: eta(nw,nlev)

! Internal
real, dimension(nlev) :: zg
real :: outf(11,nlev)              ! iri output (11,jmax)
real :: rz12,stl
real :: d(8), t(2), sw(25), oarr(30)
integer :: l,j,jj,ijf,jmag,iday,mmdd,i,ii,n,k,ix,iy
integer :: instance,iostatus
integer :: iday_iri
character(len=1024) :: iri90_dir
logical :: jf(12)
data sw/25*1./

!
! Set data paths
!
data_dir=data_path
! 
! Define input variables to be used by GLOW
!
glong = glong_in
glat = glat_in
sza = sza_in

f107 = solar_geomags(1)
f107a = solar_geomags(2)
f107p = solar_geomags(3)

iscale = switches(1)
xuvfac = switches(2)
kchem = switches(3)
jlocal = switches(4)

itail = auroral(1)
fmono = auroral(2)
emono = auroral(3)
ef = auroral(4)
ec = auroral(5)
!
! Call CGLOW_INIT (module CGLOW) to set array dimensions and allocate use-associated variables:
! (This was formerly done using common blocks, including common block /cglow/.)
!
jmax=nlev
if (.NOT. allocated(zz)) then
    call cglow_init
end if
!
! Call EGRID to set up electron energy grid:
!
call egrid (ener, del, nbins)
!
! Define background
!
zg = bks(:,1)
zz = zg*1e5 !km to cm
zo = bks(:,2)
zo2 = bks(:,3)
zn2 = bks(:,4)
ztn = bks(:,5)
ze = bks(:,6)
zte = bks(:,7)
zti = bks(:,8)
zno = bks(:,9)
znd = bks(:,10)
zns = bks(:,11)
zxden(3,:) = bks(:,12)
zxden(6,:) = bks(:,13)
zxden(7,:) = bks(:,14)
!
! Auroral inputs?
!
phitop(:) = 0.
if (ef>0.001 .and. ec>1.) call maxt(ef,ec,ener,del,nbins,itail,fmono,emono,phitop)
!
! Call GLOW
!
call glow
if (ierr > 0) write(6,"('glow_drv: IERR = ',i1)") ierr
if (abs(efrac) > 0.2) write(6,"('glow_drv:  EFRAC =',f5.2,'  SZA =',f5.2,'  DIP=',f5.2)") efrac,sza,dip
!
! output arrays
!
do ix=1,nw
    do iy=1,nlev
      if (isnan(zeta(ix,iy)) .or. zeta(ix,iy) < 0) then
        zeta(ix,iy) = 0. 
      endif
    enddo
enddo
do ix=1,nc
    do iy=1,nlev
      if (isnan(zlbh(ix,iy)) .or. zlbh(ix,iy) < 0) then
        zlbh(ix,iy) = 0.
      endif
      if (isnan(zceta(ix,5,iy)) .or. zceta(ix,5,iy) < 0) then
        zceta(ix,5,iy) = 0.
      endif
    enddo
enddo
!
! Outputs
!
!xden = zxden
eta = zeta
ceta = zceta
!
! Deallocate all variables
!
call cglow_final
  
end subroutine glow1d

end module