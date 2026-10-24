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


subroutine bands

  use cglow, only: jmax, nc, aglw, zeta, zlbh

  implicit none

  integer :: k,m
  real :: vp_pop(nc)
  data vp_pop /0.068,0.132,0.168,0.18,0.158,0.105,0.07,0.0,0.0,0.0/ ! Determined from GOLD (Ajello et al. 2020)

  do k=1,jmax
    do m=1,nc
      zlbh(m,k)=aglw(4,3,k)*vp_pop(m)
    end do
  end do

  return

end subroutine bands
