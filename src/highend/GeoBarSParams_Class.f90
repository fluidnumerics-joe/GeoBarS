! GeoBarSParams_Class.f90
! 
! Copyright 2018 Joseph Schoonover <joeschoonover@fluidnumerics.com>, Fluid Numerics, LLC
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
 
 
MODULE GeoBarSParams_Class

 USE ModelPrecision
 USE CommonRoutines
 USE ConstantsDictionary
 

 IMPLICIT NONE


    TYPE GeoBarSParams
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       INTEGER       :: mInnerIters
       INTEGER       :: pcIterates
       REAL(prec)    :: tolerance
       REAL(prec)    :: pcTolerance
       INTEGER       :: pcDegree

       ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       CHARACTER(50) :: PeaceMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: nXElem
       INTEGER       :: nYElem
       INTEGER       :: nPlot
       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale

       ! GeoBarSParameters
       REAL(prec)    :: cDrag
       REAL(prec)    :: f0
       REAL(prec)    :: betaX
       REAL(prec)    :: betaY 
      


       CONTAINS

       PROCEDURE :: Build => BuildParams_GeoBarS

    END TYPE GeoBarSParams 
 

 CONTAINS


 SUBROUTINE BuildParams_GeoBarS( thisParam )
 ! S/R BuildParams
 !
 !
 ! ==============================================!
 ! DECLARATIONS
   CLASS( GeoBarSParams ), intent(out) :: thisParam
   ! LOCAL
   INTEGER :: nUnit
   
       ! SolverCriteria
       INTEGER       :: MaximumIterates
       INTEGER       :: mInnerIters
       INTEGER       :: pcIterates
       REAL(prec)    :: tolerance
       REAL(prec)    :: pcTolerance
       INTEGER       :: pcDegree

       ! SpaceManagement
       CHARACTER(50) :: SpecMeshFile
       CHARACTER(50) :: PeaceMeshFile
       INTEGER       :: polyDeg
       INTEGER       :: nXElem
       INTEGER       :: nYElem
       INTEGER       :: nPlot
!       REAL(prec)    :: dxPlot
       REAL(prec)    :: xScale
       REAL(prec)    :: yScale

       ! GeoBarSParameters
       REAL(prec)    :: cDrag
       REAL(prec)    :: f0
       REAL(prec)    :: betaX
       REAL(prec)    :: betaY

      ! Set the namelists
      NAMELIST / SolverCriteria / MaximumIterates, mInnerIters, pcIterates, tolerance, pcTolerance, pcDegree

      NAMELIST / SpaceManagement / SpecMeshFile, PeaceMeshFile, polyDeg, nXElem, nYElem, nPlot, xScale, yScale
 
      NAMELIST / GeoBarSParameters / cDrag, f0, betaX, betaY

      ! Set the default parameters
      ! SolverCriteria
      MaximumIterates = 500   ! Max number of conjugate gradient iterations
      mInnerIters     = 50
      pcIterates      = 50     !
      tolerance       = 10.0_prec**(-8) ! conjugate gradient residual tolerance
      pcTolerance     = 10.0_prec**(-4) 
      pcDegree        = 1
      ! SpaceManagement
      SpecMeshFile = nada 
      PeaceMeshFile = nada
      polyDeg      = 5
      nXElem       = 5
      nYElem       = 5
      nPlot        = 10
      xScale       = 1.0_prec*10.0_prec**(6)
      yScale       = 3.0_prec*10.0_prec**(6)
      ! GeoBarSParameters
      cDrag = 0.01_prec
      f0    = 1.0_prec*10.0_prec**(-4)
      betaX = ZERO
      betaY = 1.0_prec*10.0_prec**(-11)

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params')
         READ( UNIT = nUnit, NML = SolverCriteria )
         READ( UNIT = nUnit, NML = SpaceManagement )
         READ( Unit = nUnit, NML = GeoBarSParameters )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SolverCriteria )
      WRITE( UNIT = *, NML = SpaceManagement )
      WRITE( UNIT = *, NML = GeoBarSParameters )

      ! Fill in the data structure
      ! SolverCriteria
      thisParam % MaximumIterates = MaximumIterates
      thisParam % mInnerIters     = mInnerIters
      thisParam % pcIterates      = pcIterates
      thisParam % tolerance       = tolerance
      thisParam % pcTolerance     = pcTolerance
      thisParam % pcDegree        = 1
      
      ! SpaceManagement
      thisParam % SpecMeshFile = SpecMeshFile
      thisParam % PeaceMeshFile = PeaceMeshFile
      thisParam % polyDeg = polyDeg
      thisParam % nXElem  = nXElem
      thisParam % nYElem  = nYElem
      thisParam % nPlot   = nPlot
      thisParam % dxPlot  = 2.0_prec/REAL(nPlot,prec)
      thisParam % xScale  = xScale
      thisParam % yScale  = yScale

      ! GeoBarSParameters
      thisParam % cDrag = cDrag
      thisParam % f0    = f0
      thisParam % betaX = betaX
      thisParam % betaY = betaY


 END SUBROUTINE BuildParams_GeoBarS

END MODULE GeoBarSParams_Class
