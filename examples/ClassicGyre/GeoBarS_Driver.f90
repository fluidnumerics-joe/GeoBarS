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
   CALL SetupBoundaryConditions( ellipticSolver )
   CALL SetEkmanPumping( ellipticSolver )

   CALL ellipticSolver % Solve( ioerr )


   CALL ellipticSolver % WritePickup( )
   CALL ellipticSolver % WriteTecplot( )
   CALL ellipticSolver % WriteResidual( )

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
 SUBROUTINE SetupBoundaryConditions( myCGSEM )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   TYPE( GeoBarS ), INTENT(inout)    :: myCGSEM
   ! LOCAL
   INTEGER    :: i, j, k, e1, e2, sID
   REAL(prec) :: x, y


      DO k = 1, myCGSEM % mesh % nEdges


         e2 = myCGSEM % mesh % edges(k)  % elementIDs(2)
         IF( e2 == DIRICHLET_INFLOW )THEN
         
            e1 = myCGSEM % mesh % edges(k) % elementIDs(1)
            sID = myCGSEM % mesh % edges(k) % elementSides(1)
            IF( sID == East .OR. sID == West )THEN

               i = myCGSEM % mesh % sideMap(sID) 
               DO j = 0, myCGSEM % N!
                   x = myCGSEM % mesh % geom(e1) % x(i,j)
                   y = myCGSEM % mesh % geom(e1) % y(i,j)
                   myCGSEM % psi(i,j,e1) = InflowStreamFunction( x, y ) 
               ENDDO 

            ELSE ! south or north sides

               j = myCGSEM % mesh % sideMap(sID)
               DO i = 0,myCGSEM % N 
                  x = myCGSEM % mesh % geom(e1) % x(i,j)
                  y = myCGSEM % mesh % geom(e1) % y(i,j)
                  myCGSEM % psi(i,j,e1) = InflowStreamFunction( x, y ) 
               ENDDO

            ENDIF
            
            myCGSEM % mesh % edges(k) % elementIDs(2) = DIRICHLET
         ENDIF
         
           
         IF( myCGSEM % mesh % edges(k) % elementIDs(2) == DIRICHLET_INFLOW )THEN
            PRINT*, 'Bad Logic'
         ENDIF

      ENDDO
      

 END SUBROUTINE SetupBoundaryConditions
!
!
!
 FUNCTION InflowStreamFunction( x, y ) RESULT( PSI )
 !
 ! =============================================================================================== !
   IMPLICIT NONE  
   REAL(prec) :: x, y
   REAL(prec) :: PSI 

      PSI = ZERO

 END FUNCTION InflowStreamFunction
!
!
!

END PROGRAM GeoBarS_Driver
