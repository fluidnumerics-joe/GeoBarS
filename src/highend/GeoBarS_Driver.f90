! GeoBarS_Driver.f90
! 
! Copyright 2018 Joseph Schoonover <joeschoonover@fluidnumerics.com>, Fluid Numerics, LLC
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

PROGRAM GeoBarS_Driver

USE ModelPrecision

USE QuadMesh_Class

USE GeoBarS_Class


IMPLICIT NONE

   TYPE( GeoBarS )       :: ellipticSolver
   INTEGER :: ioerr

   CALL ellipticSolver % Build(  )

   CALL SetCoriolisParameter( ellipticSolver )
   CALL SetBathymetry( ellipticSolver )
!   CALL SetBoundaryConditions( ellipticSolver )
   CALL SetEkmanPumping( ellipticSolver )

   CALL ellipticSolver % Solve( ioerr )


   CALL ellipticSolver % WritePickup( )
   CALL ellipticSolver % WriteTecplot( )

   CALL ellipticSolver % Trash( )

CONTAINS

 SUBROUTINE SetCoriolisParameter( myCGSEM )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( GeoBarS ), INTENT(inout)    :: myCGSEM
   ! LOCAL
   INTEGER    :: i, j, k
   REAL(prec) :: x, y, f0, beta

      
      f0  = myCGSEM % params % f0
      beta = myCGSEM % params % betaY

      DO k = 1, myCGSEM % mesh % nElems
         DO j = 0, myCGSEM % N
            DO i = 0, myCGSEM % N

               x = myCGSEM % mesh % geom(k) % x(i,j)
               y = myCGSEM % mesh % geom(k) % y(i,j)
               myCGSEM % f(i,j,k) = f0 + beta*y

            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE SetCoriolisParameter
!
!
!
  SUBROUTINE SetBathymetry( myCGSEM )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( GeoBarS ), INTENT(inout)    :: myCGSEM
   ! LOCAL
   INTEGER    :: i, j, k
   REAL(prec) :: x, y


      DO k = 1, myCGSEM % mesh % nElems
         DO j = 0, myCGSEM % N
            DO i = 0, myCGSEM % N
               
               x = myCGSEM % mesh % geom(k) % x(i,j)
               y = myCGSEM % mesh % geom(k) % y(i,j)
               myCGSEM % h(i,j,k) = 1000.0_prec

            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE SetBathymetry
!
!
!
 SUBROUTINE SetEkmanPumping( myCGSEM )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( GeoBarS ), INTENT(inout)    :: myCGSEM
   ! LOCAL
   INTEGER    :: i, j, k
   REAL(prec) :: x, y, f, Ly

      Ly  = myCGSEM % params % yScale

      DO k = 1, myCGSEM % mesh % nElems
         DO j = 0, myCGSEM % N
            DO i = 0, myCGSEM % N

               x = myCGSEM % mesh % geom(k) % x(i,j)
               y = myCGSEM % mesh % geom(k) % y(i,j)
               f = myCGSEM % f(i,j,k)

               myCGSEM % source(i,j,k) = -f*0.001_prec*( sin(pi*y/Ly)  )

            ENDDO
         ENDDO
      ENDDO

 END SUBROUTINE SetEkmanPumping
!
!
! 
 FUNCTION InflowStreamFunction( x, y )
 !
 ! =============================================================================================== !
   IMPLICIT NONE  
   REAL(prec) :: x, y
   REAL(prec) :: InflowStreamFunction 

      InflowStreamFunction = ZERO

 END FUNCTION InflowStreamFunction
!
!
!

END PROGRAM GeoBarS_Driver
