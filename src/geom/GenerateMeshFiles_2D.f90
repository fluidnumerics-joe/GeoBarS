! GenerateMeshFiles_2D.f90
! 
! Copyright 2018 Joseph Schoonover <joeschoonover@fluidnumerics.com>, Fluid Numerics, LLC
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


PROGRAM GenerateMeshFiles_2D

! src/common/
USE ModelPrecision
! src/nodal/
USE NodalStorage_Class
! src/geom/
USE QuadMesh_Class
! src/boundary/
USE BoundaryCommunicator_Class
! src/highend/shallowwater
USE Params_Class

 IMPLICIT NONE
 
   TYPE( RunParams )            :: params
   TYPE( QuadMesh )             :: mesh
   TYPE( BoundaryCommunicator ) :: extComm
   TYPE( NodalStorage )         :: dgStorage
   INTEGER                      :: N, nPlot, nBe, iEdge, e2, quadType
   
      CALL params % Build( )
      N     = params % polyDeg
      nPlot = params % nPlot
      quadType = params % QuadType
      
      CALL dGStorage % Build( N, nPlot, quadType, DG )
      
      IF( TRIM( params % SpecMeshFile ) == nada )THEN
         PRINT*,' Loading default mesh.'
         CALL mesh % LoadDefaultMesh( dgStorage % interp, &
                                      params % nXelem, &
                                      params % nYelem )
      ELSE
      ! Builds the lateral mesh
         PRINT*, 'Reading mesh from '//trim(params % SpecMeshFile)//'.'
         CALL  mesh % ReadSpecMeshFile( dgStorage % interp, &
                                        params % SpecMeshFile )
      ENDIF
      
      nBe = 0
      DO iEdge = 1, mesh % nEdges
         e2 = mesh % edges(iEdge) % ElementIDs(2)
         IF( e2 <= 0 )THEN
            nBe = nBe + 1
         ENDIF
      ENDDO
      CALL extComm % Initialize( nBe )
      
      nBe = 0
      DO iEdge = 1, mesh % nEdges
         e2 = mesh % edges(iEdge) % ElementIDs(2)
         IF( e2 <= 0 )THEN
            nBe = nBe + 1
            extComm % boundaryIDs(nBe) = iEdge
            extComm % extElemIDs(nBe)  = e2
            mesh % edges(iEdge) % boundaryID = nBe
         ENDIF

      ENDDO

      CALL mesh % WriteTecplot( )
      CALL mesh % WritePeaceMeshFile( params % PeaceMeshFile )
      CALL extComm % WritePickup( 'ExtComm' )
    
      CALL mesh % Trash( )
      CALL extComm % Trash( )

 
END PROGRAM GenerateMeshFiles_2D
