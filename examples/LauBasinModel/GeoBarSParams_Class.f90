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
       
       ! Topography
       REAL(prec) :: xS1, xS2, xN1, xN2, xC1, xC2, xC3
       REAL(prec) :: yS, yN, yC1, yC2, yC3, xI, u0
       REAL(prec) :: H, h1, h2, l1, l2, dh1, dh2, dh3
       REAL(prec) :: lx1, lx2, lx3, ly1, ly2, ly3
       REAL(prec) :: yshift, xshift, southernCutOff


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
       
       ! Topography
       REAL(prec) :: xS1, xS2, xN1, xN2, xC1, xC2, xC3
       REAL(prec) :: yS, yN, yC1, yC2, yC3, xI, u0
       REAL(prec) :: H, h1, h2, l1, l2, dh1, dh2, dh3
       REAL(prec) :: lx1, lx2, lx3, ly1, ly2, ly3
       REAL(prec) :: yshift, xshift, southernCutOff

      ! Set the namelists
      NAMELIST / SolverCriteria / MaximumIterates, mInnerIters, pcIterates, tolerance, pcTolerance, pcDegree

      NAMELIST / SpaceManagement / SpecMeshFile, PeaceMeshFile, polyDeg, nXElem, nYElem, nPlot, xScale, yScale
 
      NAMELIST / GeoBarSParameters / cDrag, f0, betaX, betaY
      
      NAMELIST / Topography / xS1, xS2, xN1, xN2, xC1, xC2, xC3, yS, yN, yC1, yC2, yC3, xI, u0, & 
                              H, h1, h2, l1, l2, dh1, dh2, dh3, lx1, lx2, lx3, ly1, ly2, ly3, &
                              yshift, xshift, southernCutOff

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
      ! Topography
      xs1 = HALF*xScale - 10.0_prec**(5)
      xs2 = xs1 + 2.0_prec*10.0_prec**(5)
      xn1 = xs1
      xn2 = xn1 + 6.0_prec*10.0_prec**(5)
      yn = 0.75_prec*yScale
      ys = yn - 15.0_prec*10.0_prec**(5)
      H = 3500.0_prec
      h1 = 700.0_prec
      h2 = 700.0_prec
      l1 = 5.0_prec*10.0_prec**(4)
      l2 = 5.0_prec*10.0_prec**(4)
      !  > Northern Passage
      yc1 = yn - 3.5_prec*10.0_prec**(5)
      xc1 = xs1 + 6.9_prec*10.0_prec**(5)
      lx1 = 5.0_prec*10.0_prec**(4)
      ly1 = 5.0_prec*10.0_prec**(4)
      dh1 = 1300.0_prec
      !  > Southern Passage
      yc2 = yn - 8.5_prec*10.0_prec**(5)
      xc2 = xs1 + 1.6_prec*10.0_prec**(5)
      lx2 = 5.0_prec*10.0_prec**(4)
      ly2 = 5.0_prec*10.0_prec**(4)
      dh2 = 800.0_prec
      !  > Southern Outlet
      yc3 = yn - 13.5_prec*10.0_prec**(5)
      xc3 = xs1
      lx3 = 4.0_prec*10.0_prec**(4)
      ly3 = 5.0_prec*10.0_prec**(4)
      dh3 = 900.0_prec
      !
      xI = 1400000.0_prec
      u0 = 0.1_prec
      yshift = 0.0_prec
      xshift = 0.0_prec
      southernCutOff = 0.0_prec

      ! Reading in the namelist FILE

      OPEN( UNIT = NEWUNIT(nUnit), FILE = 'runtime.params')
         READ( UNIT = nUnit, NML = SolverCriteria )
         READ( UNIT = nUnit, NML = SpaceManagement )
         READ( Unit = nUnit, NML = GeoBarSParameters )
         READ( Unit = nUnit, NML = Topography )
      CLOSE( UNIT = nUnit ) 

      ! Sanity check - PRINT the results of the namelist READ to screen
      WRITE( UNIT = *, NML = SolverCriteria )
      WRITE( UNIT = *, NML = SpaceManagement )
      WRITE( UNIT = *, NML = GeoBarSParameters )
      WRITE( UNIT = *, NML = Topography )

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
      
      ! Topography
      thisParam % xs1 = xs1
      thisParam % xs2 = xs2
      thisParam % xn1 = xn1
      thisParam % xn2 = xn2
      thisParam % yn  = yn
      thisParam % ys  = ys
      thisParam % H   = H
      thisParam % h1  = h1
      thisParam % h2  = h2
      thisParam % l1  = l1
      thisParam % l2  = l2
      thisParam % yshift = yshift
      thisParam % xshift = xshift
      thisParam % southernCutOff = southernCutOff
      !  > Northern Passage
      thisParam % yc1 = yc1
      thisParam % xc1 = xc1
      thisParam % lx1 = lx1
      thisParam % ly1 = ly1
      thisParam % dh1 = dh1
      !  > Southern Passage
      thisParam % yc2 = yc2
      thisParam % xc2 = xc2
      thisParam % lx2 = lx2
      thisParam % ly2 = ly2
      thisParam % dh2 = dh2
      !  > Southern Outlet
      thisParam % yc3 = yc3
      thisParam % xc3 = xc3
      thisParam % lx3 = lx3
      thisParam % ly3 = ly3
      thisParam % dh3 = dh3
      thisParam % xI  = xI
      thisParam % u0  = u0

 END SUBROUTINE BuildParams_GeoBarS

END MODULE GeoBarSParams_Class
