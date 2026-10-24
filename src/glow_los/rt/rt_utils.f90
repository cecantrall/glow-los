module rt_utils

implicit none

contains

subroutine read_o2fuvxsect(wavelength,o2xsect,o2fuvxsect_file)
  use netcdf
  
  ! This is the name of the data file to be read
  character (len = 200), intent(in) :: o2fuvxsect_file 
  
  ! Reading 1D data for two parameters
  integer, parameter :: bin_size = 1200
  real, intent(out) :: wavelength(bin_size), o2xsect(bin_size)
  
  ! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid, istat

  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file.
  istat = nf90_open(trim(o2fuvxsect_file)//'o2absxsectfuv.nc', NF90_NOWRITE, ncid)
  
  ! Get the varid of the wavelength variable, based on its name.
  istat = nf90_inq_varid(ncid, "WAVELENGTH", varid) 
  ! Read the data.
  istat = nf90_get_var(ncid, varid, wavelength)
  
  ! Get the varid of the cross section variable, based on its name.
  istat = nf90_inq_varid(ncid, "O2XSECT", varid)
  ! Read the data.
  istat = nf90_get_var(ncid, varid, o2xsect)
  
  ! Close the file, freeing all resources.
  istat = nf90_close(ncid) 
  
end subroutine read_o2fuvxsect

! -------------------------------------------------------------------------------------------

subroutine interpolate_value( xData, yData, xVal, yVal )
! Inputs: xData = a vector of the x-values of the data to be interpolated
!         yData = a vector of the y-values of the data to be interpolated
!         xVal  = a vector of the x-values where interpolation should be performed
! Output: yVal  = a vector of the resulting interpolated values
 
! Define inputs and outputs
real, intent(in) :: xData(:), yData(:), xVal
real, intent(out) :: yVal
! Define internal variables
integer :: dataIndex
real :: weight
 
dataIndex = 1
do while (xData(dataIndex) < xVal)
   dataIndex = dataIndex+1
enddo
weight = (xVal - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex));
yVal = (1.0-weight)*yData(dataIndex) + weight*yData(dataIndex+1);
 
end subroutine interpolate_value

! -------------------------------------------------------------------------------------------

subroutine interpolate_1D( xData, yData, xVal, yVal )
  ! Inputs: xData = a vector of the x-values of the data to be interpolated
  !         yData = a vector of the y-values of the data to be interpolated
  !         xVal  = a vector of the x-values where interpolation should be performed
  ! Output: yVal  = a vector of the resulting interpolated values
  
  ! Define inputs and outputs
  real, intent(in) :: xData(:), yData(:), xVal(:)
  real, intent(out) :: yVal(:)
  ! Define internal variables
  integer :: inputIndex, dataIndex
  real :: weight
  
  do inputIndex = 1, size(xVal)
    dataIndex = 1
    do while (xData(dataIndex) < xVal(inputIndex))
      dataIndex = dataIndex+1
    enddo
    weight = (xVal(inputIndex) - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex));
    yVal(inputIndex) = (1.0-weight)*yData(dataIndex) + weight*yData(dataIndex+1);
  end do
  
end subroutine interpolate_1D

! ---------------------------------------------------------------------------------------------

subroutine o2_scd(scdo2, o2, dpath, nlev)

    ! Define I/O
    integer, intent(in) :: nlev
    real, intent(in) :: dpath(nlev) ! cm
    real, intent(in) :: o2(nlev) ! cm
    real, intent(out) :: scdo2(nlev) ! cm^-2

    ! looping
    integer :: k
    
    do k = 1,nlev
      scdo2(k) = dot_product(o2(k:nlev),dpath(k:nlev))
    end do

end subroutine o2_scd

! ---------------------------------------------------------------------------------------------

subroutine o_scd(scdo, o, dpath, nlev)

! Define I/O
    integer, intent(in) :: nlev
    real, intent(in) :: dpath(nlev), o(nlev)
    real, intent(out) :: scdo(nlev)
! looping
integer :: k
 
do k = 1,nlev
    scdo(k) = dot_product(o(k:nlev),dpath(k:nlev))
end do
 
end subroutine o_scd

! --------------------------------------------------------------------------------------------

subroutine oxsect(xox0,xox1356,xox1359)
 
! Define variables for atomic oxygen cross sections
real, intent(out) :: xox1356, xox1359
real, intent(in) :: xox0
 
xox1356=xox0 * 0.75 / 2
xox1359=xox0 * 0.15 / 2
 
end subroutine oxsect

! ---------------------------------------------------------------------------------------------
! Below are integral routines taken from redister that have not yet been included in line of sight integration. These need to be eventually tested and included

!SUBROUTINE INTEGE(DX,Y,Z,N)
!      IMPLICIT REAL*8 (A-H,O-Z)
!********1*********2*********3*********4*********5*********6*********7**
! EVALUATE Z(I) = INTEGRAL FROM X(1) TO X(I) OF Y, ASSUMING THAT Y
! VARIES EXPONENTIALLY WITH X.
!********1*********2*********3*********4*********5*********6*********7**
!      DIMENSION DX(400),Y(400),Z(400)
!      Z(1) = 0.0D0
!      DO I = 2,N
!         DX = X(I) - X(I - 1)
!         EA = Y(I - 1)
!         EB = Y(I)
!         IF(EA .NE. EB .AND. EA .GT. 1.0D-20 .AND.
!     &      EB .GT. 1.0D-20) THEN
!           Z(I) = Z(I - 1) + DX * (EB - EA) / LOG(EB / EA)
!         ELSE
!           Z(I) = Z(I - 1) + DX * (EB + EA) / 2.0D0
!         END IF
!      END DO
!      RETURN
!END SUBROUTINE INTEGE

!SUBROUTINE INTEGS(X,Y,Z,N,RESULT)
!      IMPLICIT REAL*8 (A-H,O-Z)
!********1*********2*********3*********4*********5*********6*********7**
! EVALUATE RESULT = INTEGRAL FROM X(1) TO X(N) OF Y * Z, ASSUMING THAT Y
! VARIES EXPONENTIALLY WITH X, WHILE Z IS LINEAR IN X.
!********1*********2*********3*********4*********5*********6*********7**
!      DIMENSION X(400),Y(400),Z(400)
!      RESULT = 0.0D0
!      DO I = 2,N
!         DX = X(I) - X(I - 1)
!         DY = Y(I) - Y(I - 1)
!         ZA = Z(I - 1)
!         ZB = Z(I)
!         EA = EXP(-Y(I - 1))
!         EB = EXP(-Y(I))
!         IF(DY .GE. 0.1D0 .AND. ZA .GT. 1.0D-20 .AND.
!     &      ZB .GT. 1.0D-20) THEN
!           RESULT = RESULT + DX * (ZA * EA - ZB * EB - (ZA - ZB) *
!     &              (EA - EB) / DY) / DY
!         ELSE
!           RESULT = RESULT + DX * (ZA * EA + ZB * EB) / 2.0D0
!         END IF
!      END DO
!      RETURN
!      END

!SUBROUTINE INTEG(X,Y,N,RESULT)
!      IMPLICIT REAL*8 (A-H,O-Z)
!********1*********2*********3*********4*********5*********6*********7**
! EVALUATE RESULT = INTEGRAL FROM X(1) TO X(N) OF Y,  USING A CUBIC
! SPLINE FIT TO Y, PLUS THE SPLINE QUADRATURE FORMULA OF FORSYTHE,
! MALCOLM, AND MOLER, P. 90.
!********1*********2*********3*********4*********5*********6*********7**
!      DIMENSION X(400),Y(400),B(400),C(400),D(400)
!      CALL SPLINE(N,X,Y,B,C,D)
!      RESULT = 0.0D0
!      DO I = 2,N
!         DX = X(I) - X(I - 1)
!         RESULT = RESULT + DX * ((Y(I) + Y(I - 1)) / 2.0D0 -
!     &                     DX * DX * (C(I) + C(I - 1)) / 12.0D0)
!      END DO
!      RETURN
!      END

end module rt_utils