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
   REAL(prec) :: x, y, x0


      DO k = 1, myCGSEM % mesh % nElems
         DO j = 0, myCGSEM % N
            DO i = 0, myCGSEM % N
               
               x = myCGSEM % mesh % geom(k) % x(i,j) - myCGSEM % params % xshift
               y = myCGSEM % mesh % geom(k) % y(i,j) - myCGSEM % params % yshift

               myCGSEM % h(i,j,k) = myCGSEM % params % H - &
                                    myCGSEM % params % dh1*&
                                    exp( -(x-myCGSEM % params % xC1)**2/(2.0_prec*myCGSEM % params % lx1**2) )*&
                                    0.5_prec*(1.0_prec -tanh( (y-myCGSEM % params % yC3)/(myCGSEM % params % ly3 ) ))*&
                                    ( 1.0_prec - exp( &
                                                 -(y-myCGSEM % params % yc1)**2/(2.0_prec*myCGSEM % params % ly1**2) ) )

               myCGSEM % psi(i,j,k)  = EastBoundaryCondition( x, y, 0.0_prec, myCGSEM % params % u0 ) 


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

               myCGSEM % source(i,j,k) = ZERO !1.0D-17*( sin(pi*y/Ly)  )* &
                    !0.5_prec*(1.0_prec + tanh((x - 3.5D6)/(0.5D6) )   )

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
   REAL(prec) :: x, y, Ly

      Ly  = myCGSEM % params % yScale
      
      DO k = 1, myCGSEM % mesh % nEdges


         e2 = myCGSEM % mesh % edges(k)  % elementIDs(2)
         IF( e2 == DIRICHLET )THEN
         
            e1 = myCGSEM % mesh % edges(k) % elementIDs(1)
            sID = myCGSEM % mesh % edges(k) % elementSides(1)
            IF( sID == East .OR. sID == West )THEN

               i = myCGSEM % mesh % sideMap(sID) 
               DO j = 0, myCGSEM % N!
                   x = myCGSEM % mesh % geom(e1) % x(i,j)
                   y = myCGSEM % mesh % geom(e1) % y(i,j)
                   myCGSEM % psi(i,j,e1) = EastBoundaryCondition( x, y, myCGSEM % params % xI, &
                                                                        myCGSEM % params % u0 ) 
               ENDDO 

            ELSE ! south or north sides

               j = myCGSEM % mesh % sideMap(sID)
               DO i = 0,myCGSEM % N 
                  x = myCGSEM % mesh % geom(e1) % x(i,j)
                  y = myCGSEM % mesh % geom(e1) % y(i,j)
                   myCGSEM % psi(i,j,e1) = EastBoundaryCondition( x, y, myCGSEM % params % xI, &
                                                                        myCGSEM % params % u0 ) 
               ENDDO

            ENDIF
            
            myCGSEM % mesh % edges(k) % elementIDs(2) = DIRICHLET
            
         ENDIF

      ENDDO
      

 END SUBROUTINE SetupBoundaryConditions
!
!
!
 FUNCTION EastBoundaryCondition( x, y, xc, u0 ) RESULT( PSI )
 !
 ! =============================================================================================== !
   IMPLICIT NONE  
   REAL(prec) :: x, y, xc, u0
   REAL(prec) :: PSI 
   REAL(prec) :: L
   REAL(prec) :: y1, y2

      L  = 75000.0_prec
      y1 = 2.75_prec*10.0_prec**6
      y2 = 1.75_prec*10.0_prec**6
      
      PSI = -u0*L*( tanh( ( y - y1)/L ) + tanh( (y-y2)/L ) +  tanh( y1/L ) + tanh( y2/L )  )

 END FUNCTION EastBoundaryCondition
!
!
!
 FUNCTION NorthBoundaryCondition( x, y, xc, u0 ) RESULT( PSI )
 !
 ! =============================================================================================== !
   IMPLICIT NONE  
   REAL(prec) :: x, y, xc, u0
   REAL(prec) :: PSI 
   REAL(prec) :: L

      L  = 75000.0_prec
      PSI = EastBoundaryCondition( x, y, xc, u0 )

 END FUNCTION NorthBoundaryCondition
!
!
!
 FUNCTION WestBoundaryCondition( x, y, xc, u0 ) RESULT( PSI )
 !
 ! =============================================================================================== !
   IMPLICIT NONE  
   REAL(prec) :: x, y, xc, u0
   REAL(prec) :: PSI 

      
      PSI = EastBoundaryCondition( x, y, xc, u0 )

 END FUNCTION WestBoundaryCondition
!
!
!
 FUNCTION SouthBoundaryCondition( x, y, xc, u0 ) RESULT( PSI )
 !
 ! =============================================================================================== !
   IMPLICIT NONE  
   REAL(prec) :: x, y, xc, u0
   REAL(prec) :: PSI 
   REAL(prec) :: L

      L  = 75000.0_prec
      
      PSI = EastBoundaryCondition( x, y, xc, u0 )

 END FUNCTION SouthBoundaryCondition
!
!
!

END PROGRAM GeoBarS_Driver
