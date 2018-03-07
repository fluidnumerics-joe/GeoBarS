! GeoBarS_Class.f90
! 
! Copyright 2018 Joseph Schoonover <joeschoonover@fluidnumerics.com>, Fluid Numerics, LLC
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE GeoBarS_Class

! src/common/
USE ModelPrecision
USE ConstantsDictionary
USE CommonRoutines
! src/spectralops/
USE NodalStorage_Class
! src/geom/
USE QuadMesh_Class
USE Edge_Class
USE Node_Class
USE MappedGeometry_2D_Class
! src/highend/GeoBarS/
USE GeoBarSParams_Class


IMPLICIT NONE

      TYPE GeoBarS
         INTEGER                   :: maxIters, nPlot, N, nresi
         REAL(prec)                :: tolerance
         REAL(prec), ALLOCATABLE   :: resi(:)
         TYPE( GeoBarSParams )     :: params
         TYPE( NodalStorage )      :: cgStorage
         TYPE( QuadMesh )          :: mesh      
         REAL(prec), ALLOCATABLE   :: psi(:,:,:)
         REAL(prec), ALLOCATABLE   :: source(:,:,:)
         REAL(prec), ALLOCATABLE   :: f(:,:,:)
         REAL(prec), ALLOCATABLE   :: h(:,:,:), u(:,:,:), v(:,:,:)

         CONTAINS

         PROCEDURE :: Build                => Build_GeoBarS
         PROCEDURE :: BuildQuadMesh        => BuildQuadMesh_GeoBarS
         PROCEDURE :: Trash                => Trash_GeoBarS
       
         PROCEDURE :: Mask                          => Mask_GeoBarS
         PROCEDURE :: UnMask                        => UnMask_GeoBarS
         PROCEDURE :: GlobalSum                     => GlobalSum_GeoBarS
!         PROCEDURE :: FluxDivergence                => FluxDivergence_GeoBarS
!         PROCEDURE :: SetThisDirichletEdge          => SetThisDirichletEdge_GeoBarS
!         PROCEDURE :: SetDirichletBoundaryCondition => SetDirichletBoundaryCondition_GeoBarS
!         PROCEDURE :: CalculateGradient             => CalculateGradient_GeoBarS
         PROCEDURE :: MatrixAction                  => MatrixAction_GeoBarS
!         PROCEDURE :: BoundaryFlux                  => BoundaryFlux_GeoBarS
         PROCEDURE :: Residual                      => Residual_GeoBarS
         PROCEDURE :: DotProduct                    => DotProduct_GeoBarS
         PROCEDURE :: Solve                         => Solve_GMRES
         
         PROCEDURE :: DiagnoseVelocity => DiagnoseVelocity_GeoBarS
         PROCEDURE :: CoarseToFine  => CoarseToFine_GeoBarS
         PROCEDURE :: WriteTecplot  => WriteTecplot_GeoBarS
         PROCEDURE :: WriteResidual => WriteResidual_GeoBarS
         PROCEDURE :: WritePickup   => WritePickup_GeoBarS
         PROCEDURE :: ReadPickup    => ReadPickup_GeoBarS
         
      END TYPE GeoBarS

      

CONTAINS

 SUBROUTINE Build_GeoBarS( myCGSEM )
 ! S/R Build
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(inout) :: myCGSEM
   ! LOCAL
   INTEGER                 :: N, nPlot, nPairs

      
      CALL myCGSEM % params % Build( )
      N      = myCGSEM % params % polyDeg
      nPlot  = myCGSEM % params % nPlot

      myCGSEM % tolerance = myCGSEM % params % tolerance
      myCGSEM % maxIters  = myCGSEM % params % maximumIterates
      myCGSEM % N         = N
      myCGSEM % nPlot     = nPlot
      myCGSEM % nresi     = myCGSEM % maxIters*myCGSEM % params % mInnerIters
     
      CALL myCGSEM % cgStorage % Build( N, nPlot, GAUSS_LOBATTO, CG )
      CALL myCGSEM % BuildQuadMesh( )

      ! Allocate space for the solution, source term, and diffusivity coefficient
      ALLOCATE( myCGSEM % psi(0:N, 0:N, 1:myCGSEM % mesh % nElems), &
                myCGSEM % source(0:N, 0:N, 1:myCGSEM % mesh % nElems), &
                myCGSEM % f(0:N, 0:N, 1:myCGSEM % mesh % nElems), &
                myCGSEM % h(0:N, 0:N, 1:myCGSEM % mesh % nElems), &
                myCGSEM % u(0:N, 0:N, 1:myCGSEM % mesh % nElems), &
                myCGSEM % v(0:N, 0:N, 1:myCGSEM % mesh % nElems), &
                myCGSEM % resi(0:myCGSEM % nresi) )
                
      myCGSEM % psi    = ZERO
      myCGSEM % source = ZERO
      myCGSEM % h      = ONE
      myCGSEM % f      = ZERO
      myCGSEM % u      = ZERO
      myCGSEM % v      = ZERO
      myCGSEM % resi   = ZERO
      
      
     ! CALL myCGSEM % ReadPickup( )
      
 END SUBROUTINE Build_GeoBarS
!
!
!
 SUBROUTINE BuildQuadMesh_GeoBarS( myCGSEM )
 ! S/R BuildQuadMesh
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(inout) :: myCGSEM

      ! Builds the lateral mesh
      PRINT*, 'Reading 2-D mesh from '//trim(myCGSEM % params % SpecMeshFile)//'.'
     ! CALL myCGSEM % mesh % ReadPeaceMeshFile( myCGSEM % params % PeaceMeshFile )
      
      CALL myCGSEM % mesh % ReadSpecMeshFile( myCGSEM % cgStorage % interp, &
                                              myCGSEM % params % SpecMeshFile )
                                              
      CALL myCGSEM % mesh % ScaleTheMesh( myCGSEM % cgStorage % interp, &
                                          myCGSEM % params % xScale, &
                                          myCGSEM % params % yScale )

 END SUBROUTINE BuildQuadMesh_GeoBarS
!
!
!
 SUBROUTINE Trash_GeoBarS( myCGSEM )
 ! S/R Trash
 !
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(inout) :: myCGSEM


      ! Trash the nodal storage structure
      CALL myCGSEM % cgStorage % Trash( )

      ! Trash the geometry
      CALL myCGSEM % mesh % Trash( )
      

      ! Deallocate space for the solution, source term, and diffusivity coefficient
      DEALLOCATE( myCGSEM % psi,&
                  myCGSEM % source, &
                  myCGSEM % h, &
                  myCGSEM % u, &
                  myCGSEM % v )
                  
 END SUBROUTINE Trash_GeoBarS
!
!
!==================================================================================================!
!-------------------------------------- Type Specific ---------------------------------------------!
!==================================================================================================!
!
!
 SUBROUTINE Mask_GeoBarS( myCGSEM, u )
 ! S/R Mask_GeoBarS
 !
 !
 ! ============================================================================= !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(inout) :: myCGSEM
   REAL(prec), INTENT(inout)                :: u(0:myCGSEM % N, &
                                                 0:myCGSEM % N, &
                                                 1:myCGSEM % mesh % nElems )
   !LOCAL
   INTEGER :: N, iEdge, eID, iNode, nID, i, j, k, nEl, nEdges, nNodes
   INTEGER :: e1, e2, s1, s2, nodeType, nElForNode, nstart

      N      = myCGSEM % N
      nEl    = myCGSEM % mesh % nElems
      nEdges = myCGSEM % mesh % nEdges
      nNodes = myCGSEM % mesh % nNodes
   
      ! Mask out the solution over the edges -- excluding the corner nodes
      DO iEdge = 1, nEdges ! Loop over the edges
 
         e2 = myCGSEM % mesh % edges(iEdge) % elementIDs(2)

         IF( e2 == DIRICHLET )THEN ! check the secondary element for boundary condition flag

            ! This is a prescribed boundary condition (Dirichlet) so we mask out u along this edge   

            e1 = myCGSEM % mesh % edges(iEdge) % elementIDs(1)
            s1 = myCGSEM % mesh % edges(iEdge) % elementSides(1) 
            CALL MaskSide( e1, s1, myCGSEM % mesh, u, N, nEl )

         ELSE ! then this is not a prescribed boundary, and we choose to mask out the secondary element
         
            e2 = myCGSEM % mesh % edges(iEdge) % elementIDs(2)
            s2 = myCGSEM % mesh % edges(iEdge) % elementSides(2)
            IF( e2 > 0 )THEN ! this is an internal edge, and we should mask out the secondary element
               CALL MaskSide( e2, s2, myCGSEM % mesh, u, N, nEl )
            ENDIF

         ENDIF

      ENDDO ! iEdge, loop over the edges     

      ! At this point, the secondary internal edges and prescribed boundary condition edges have been masked out.
      ! Now we mask out all but the primary element in the corner-node lists.

      ! Mask out the corner-nodes appropriately ( Leave the primary element intact )
      DO iNode = 1, nNodes ! Loop over the corner nodes

         nElForNode = myCGSEM % mesh % nodeToElement(iNode,0,1)
         !CALL myCGSEM % mesh % nodes(iNode) % RewindElementList( )
         nodeType = myCGSEM % mesh % nodes(iNode) % nodeType

         IF( nodeType == DIRICHLET .OR. nodeType == DIRICHLET_INFLOW )THEN
            nstart = 1
         ELSE 
            nstart = 2
         ENDIF
         
         DO k = nstart, nElForNode ! while there are elements in the node-to-element list
            
            eID = myCGSEM % mesh % nodeToElement(iNode, k, 1)
            nID = myCGSEM % mesh % nodeToElement(iNode, k, 2)
            
            !CALL myCGSEM % mesh % nodes(iNode) % GetCurrentElementID( eID ) ! Get the element ID and the local node ID (1->4)
            !CALL myCGSEM % mesh % nodes(iNode) % nodeToElement % GetKey( nID ) ! Get the local node ID

            i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
            j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
 
            u(i,j,eID) = ZERO ! Mask out the solution at this element at this corner node

         ENDDO ! while we have elements in the node-to-element list

      ENDDO


 END SUBROUTINE Mask_GeoBarS
!
!
!
 SUBROUTINE MaskSide( eID, sID, mesh, u, N, nElems )
 ! S/R MaskSide
 !
 !  This subroutine masks the solutions on side "sID" of element "eID",
 !  excluding the corner nodes ( Assumes Gauss-Lobatto points are used )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)        :: eID, sID, N,nElems
   REAL(prec), INTENT(inout)  :: u(0:N,0:N,1:nElems)
   TYPE(QuadMesh), INTENT(in) :: mesh
   ! LOCAL
   INTEGER :: iS, iP, si
   
      si = ABS(sID)
      IF( si == East .OR. si == West )then ! east or west sides

         iS = mesh % sideMap(si) 
         DO iP = 1, N-1 ! Loop over the edge, excluding the corner nodes
            u(iS,iP,eID) = ZERO
         ENDDO ! iP, loop over the edge

      ELSE ! south or north sides

         iP = mesh % sideMap(si)
         DO iS = 1, N-1 ! Loop over the edge, excluding the corner nodes
            u(iS,iP,eID) = ZERO
         ENDDO

      ENDIF

 END SUBROUTINE MaskSide
!
!
!
 SUBROUTINE UnMask_GeoBarS( myCGSEM, u )
 ! S/R UnMask
 !
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(inout) :: myCGSEM
   REAL(prec), INTENT(inout)       :: u(0:myCGSEM % N, &
                                        0:myCGSEM % N, &
                                        1:myCGSEM % mesh % nElems )
   !LOCAL
   INTEGER :: N, iEdge, eID, iNode, nID, i, j, k, eID0, i0, j0
   INTEGER :: nodeType, e2, nElForNode, bID

      ! Unmask the solution over the interior edges only -- corner nodes are excluded in this step
      ! and are taken care of in the following step
      DO iEdge = 1, myCGSEM % mesh % nEdges ! Loop over the edges
 
         e2 = myCGSEM % mesh % edges(iEdge) % elementIDs(2)
         IF( e2 > 0 )THEN ! This is an interior shared edge
            CALL UnMaskSide( iEdge, myCGSEM % mesh, u, myCGSEM % N, myCGSEM % mesh % nElems ) ! unmask the solution along this edge
         ENDIF

      ENDDO ! iEdge, loop over the edges     


      ! Unmask the corner-nodes
      DO iNode = 1, myCGSEM % mesh % nNodes ! Loop over the corner nodes
         
         nElForNode = myCGSEM % mesh % nodeToElement(iNode,0,1)
         nodeType = myCGSEM % mesh % nodes(iNode) % nodeType

         IF( nodeType /= DIRICHLET .AND. nodeType /= DIRICHLET_INFLOW ) then ! this is NOT a prescribed (Dirichlet) boundary

            !CALL myCGSEM % mesh % nodes(iNode) % RewindElementList( )
            eID0 = myCGSEM % mesh % nodeToElement(iNode, 1, 1)
            nID = myCGSEM % mesh % nodeToElement(iNode, 1, 2)
            
            ! Gather the element and local node IDs of the corner which contains the unmasked solution values.
            i0 = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
            j0 = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node

            !CALL myCGSEM % mesh % nodes(iNode) % AdvanceElementList( ) ! move on to the next element which shares this node


            ! Loop over the elements in this corner-node's connectivity list
            DO k = 2, nElForNode ! while there are elements in the node-to-element list

               eID = myCGSEM % mesh % nodeToElement(iNode, k, 1)
               nID = myCGSEM % mesh % nodeToElement(iNode, k, 2)
 
               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               u(i,j,eID) = u(i0,j0,eID0) ! copy the solution from the unmasked solution

            !   CALL myCGSEM % mesh % nodes(iNode) % AdvanceElementList( )

            ENDDO ! while we have elements in the node-to-element list

         
         ENDIF

      ENDDO

 END SUBROUTINE UnMask_GeoBarS
!
!
!
 SUBROUTINE UnMaskSide( edgeID, mesh, u, N, nElems )
 ! S/R UnMaskSide
 !
 !  This subroutine unmasks the solutions on side "sID" of element "eID",
 !  excluding the corner nodes ( Assumes Gauss-Lobatto points are used )
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)        :: edgeID
   INTEGER, INTENT(in)        :: N, nElems
   REAL(prec), INTENT(inout)  :: u(0:N,0:N,1:nElems)
   TYPE(QuadMesh), INTENT(in) :: mesh
   ! LOCAL
   INTEGER :: iS, iP, eID, sID, jS, jP, inc
   REAL(prec) :: temp(0:N) ! assume nS == nP

      ! Get the primary element and side IDs
      eID = mesh % edges(edgeID) % elementIDs(1)
      sID = mesh % edges(edgeID) % elementSides(1)
 
      IF( sID == East .OR. sID == West )THEN ! east or west sides

         iS = mesh % sideMap(sID) 
         DO iP = 1, N-1 ! Loop over the edge, excluding the corner nodes
            temp(iP) = u(iS,iP,eID) ! Copy the solution in the primary element
         ENDDO ! iP, loop over the edge

      ELSE ! south or north sides

         iP = mesh % sideMap(sID)
         DO iS = 1, N-1 ! Loop over the edge, excluding the corner nodes
            temp(iS) = u(iS,iP,eID) ! Copy the solution in the primary element
         ENDDO

      ENDIF

       ! Get the secondary element and side IDs
      eID = mesh % edges(edgeID) % elementIDs(2)
      sID = mesh % edges(edgeID) % elementSides(2)
      sID = abs(sID)
      
      IF( sID == East .OR. sID == West )THEN ! east or west sides

         iS = mesh % sideMap(sID) 
         
         jP  = mesh % edges(edgeID) % start
         inc = mesh % edges(edgeID) % inc
         
         DO iP = 1, N-1 ! Loop over the edge, excluding the corner nodes
            u(iS,jP,eID) = temp(iP) ! Copy the solution to the secondary element
            jP = jP + inc ! increment jP
         ENDDO ! iP, loop over the edge

      ELSE ! south or north sides

         iP = mesh % sideMap(sID)
         jS  = mesh % edges(edgeID) % start
         inc = mesh % edges(edgeID) % inc

         DO iS = 1, N-1 ! Loop over the edge, excluding the corner nodes

            u(jS,iP,eID) = temp(iS) ! Copy the solution in the primary element
            jS = jS + inc ! increment jS

         ENDDO

      ENDIF

 END SUBROUTINE UnMaskSide
!
!
!
 SUBROUTINE GlobalSum_GeoBarS( myCGSEM, u )
 ! S/R GlobalSum
 !
 !  Adds together the shared edge and node contributions for the CGSEM method.
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(inout) :: myCGSEM
   REAL(prec), INTENT(inout)                :: u(0:myCGSEM % N, &
                                                 0:myCGSEM % N, &
                                                 1:myCGSEM % mesh % nElems )
   !LOCAL
   INTEGER    :: iEdge, eID, iNode, nID, i, j, k, nEdges, nNodes, N, nEl
   INTEGER    :: e2, nodeType, nElForNode
   REAL(prec) :: theSum

      N      = myCGSEM % N
      nEl    = myCGSEM % mesh % nElems
      nEdges = myCGSEM % mesh % nEdges
      nNodes = myCGSEM % mesh % nNodes

      DO iEdge = 1, nEdges ! loop over the edges

         e2 = myCGSEM % mesh % edges(iEdge) % elementIDs(2)
         IF( e2 > 0  )then ! This is an interior shared edge 

            CALL SumSide( iEdge, myCGSEM % mesh, u, N, nEl )  ! add the contributions from the elements which share this edge
                                                                   ! The edge sum excludes corner points
         ENDIF

      ENDDO ! iEdge, loop over the edges


      DO iNode = 1, nNodes ! loop over the corner nodes
      
         nElForNode = myCGSEM % mesh % nodeToElement(iNode,0,1)
         nodeType = myCGSEM % mesh % nodes(iNode) % nodeType

         IF( nodeType /= DIRICHLET ) then ! this node does NOT lie on a prescribed (Dirichlet) boundary

            ! *** First, compute the sum from the contributing elements *** !
            theSum = ZERO 

            ! Loop over the elements in this corner-node's connectivity list
            DO k = 1, nElForNode ! while there are elements in the node-to-element list

               eID = myCGSEM % mesh % nodeToElement(iNode,k,1)
               nID = myCGSEM % mesh % nodeToElement(iNode,k,2)

               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               theSum = theSum + u(i,j,eID) 

              ! CALL myCGSEM % mesh % nodes(iNode) % AdvanceElementList( )

            ENDDO ! while we have elements in the node-to-element list

            ! *** Now, copy the sum to the array "u" for each element which contributed to the sum *** !

            DO k = 1, nElForNode ! while there are elements in the node-to-element list

               eID = myCGSEM % mesh % nodeToElement(iNode,k,1)
               nID = myCGSEM % mesh % nodeToElement(iNode,k,2)
               
               i = myCGSEM % mesh % cornerMap(1, nID) ! Grab the index over the first computational direction associated with this corner node
               j = myCGSEM % mesh % cornerMap(2, nID) ! Grab the index over the second computational direction associated with this corner node
  
               u(i,j,eID) = theSum

            ENDDO ! while we have elements in the node-to-element list

         ENDIF

      ENDDO

 END SUBROUTINE GlobalSum_GeoBarS
!
!
!
 SUBROUTINE SumSide( edgeID, mesh, u, N, nElems )
 ! S/R SumSide
 !
 ! 
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   INTEGER, INTENT(in)        :: edgeID
   INTEGER, INTENT(in)        :: N, nElems
   REAL(prec), INTENT(inout)  :: u(0:N,0:N,1:nElems)
   TYPE(QUADMESH), INTENT(in) :: mesh
   ! LOCAL
   INTEGER :: iS, iP, k, inc, e(1:2), s(1:2)
   REAL(prec) :: temp(0:N,1:2) ! assume nS == nP
   REAL(prec) :: theSum

      ! Get the primary/secondary element and side IDs
      e = mesh % edges(edgeID) % elementIDs
      s = mesh % edges(edgeID) % elementSides  
      s(2) = abs(s(2))

      DO k = 1, 2
         IF( s(k) == east .OR. s(k) == west )then ! east or west sides

            iS = mesh % sideMap(s(k)) 
            DO iP = 1, N-1 ! Loop over the edge, excluding the corner nodes
               temp(iP,k) = u(iS,iP,e(k)) ! 
            ENDDO ! iP, loop over the edge

         ELSE ! south or north sides

            iP = mesh % sideMap(s(k))
            DO iS = 1, N-1 ! Loop over the edge, excluding the corner nodes
               temp(iS,k) = u(iS,iP,e(k)) ! 
            ENDDO

         ENDIF

      ENDDO

      ! Add the two contributions together
      k   = mesh % edges(edgeID) % start
      inc = mesh % edges(edgeID) % inc
     
      DO iS = 1, N-1

         theSum = temp(iS,1) + temp(k,2)
         temp(iS,1) = theSum
         temp(k,2) = theSum

         k = k + inc ! increment the secondary element array adress

      ENDDO

      ! copy the sum into the two contributing element edges

      DO k = 1, 2
   
         IF( s(k) == east .OR. s(k) == west )THEN ! east or west sides

            iS = mesh % sideMap(s(k)) 
            DO iP = 1, N-1 ! Loop over the edge, excluding the corner nodes
                u(iS,iP,e(k)) = temp(iP,k) !
            ENDDO ! iP, loop over the edge

         ELSE ! south or north sides

            iP = mesh % sideMap(s(k))
            DO iS = 1, N-1 ! Loop over the edge, excluding the corner nodes
               u(iS,iP,e(k)) = temp(iS,k) ! 
            ENDDO

         ENDIF

      ENDDO

 END SUBROUTINE SumSide
!
!
!
 SUBROUTINE MatrixAction_GeoBarS( myCGSEM, u, Au )
 ! S/R MatrixAction
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(GeoBarS), INTENT(inout) :: myCGSEM
   REAL(prec), INTENT(in)        :: u(0:myCGSEM % N, &
                                      0:myCGSEM % N, &
                                      1:myCGSEM % mesh % nElems )
   REAL(prec), INTENT(out)       :: Au(0:myCGSEM % N, &
                                       0:myCGSEM % N, &
                                       1:myCGSEM % mesh % nElems  )
   ! Local
   INTEGER :: i, j, k, iEl
   REAL(prec) :: duds, dudp, dudx, dudy, dF1ds, dF2dp, Fx, Fy
   REAL(prec) :: F1(0:myCGSEM % N, 0:myCGSEM % N)
   REAL(prec) :: F2(0:myCGSEM % N, 0:myCGSEM % N)

     
!$OMP PARALLEL


!$OMP DO PRIVATE( duds, dudp, dudx, dudy, Fx, Fy, F1, F2, dF1ds, dF2dp )
      DO iEl = 1, myCGSEM % mesh % nElems
      
      
         !------------- Contravariant Flux calculations ---------------!
         ! Calculate the gradient in the stream function 
         !CALL myCGSEM % CalculateGradient( iEl, u, dudx, dudy ) ! calculate the gradient in physical coordinates
         DO j = 0, myCGSEM % N
            DO i = 0, myCGSEM % N
            
               duds = 0.0_prec
               dudp = 0.0_prec
               
            !   duds = MATMUL( myCGSEM % cgStorage % dMatS, u )
            !   dudp = MATMUL( u, myCGSEM % cgStorage % dMatP )
               ! Gradient in computational space
               DO k = 0, myCGSEM % N
                  duds = duds + myCGSEM % cgStorage % dMatS(i,k)*u(k,j,iEl)
                  dudp = dudp + u(i,k,iEl)*myCGSEM % cgStorage % dMatP(k,j)
               ENDDO
               
               
               ! Calculate the gradient in physical space
               dudx = ( myCGSEM % mesh % geom(iEl) % dydp(i,j)*duds - &
                        myCGSEM % mesh % geom(iEl) % dyds(i,j)*dudp )/&
                        myCGSEM % mesh % geom(iEl) % J(i,j)
                        
               dudy = ( myCGSEM % mesh % geom(iEl) % dxds(i,j)*dudp - &
                        myCGSEM % mesh % geom(iEl) % dxdp(i,j)*duds )/&
                        myCGSEM % mesh % geom(iEl) % J(i,j)

               !Fx = myCGSEM % params % cDrag*dudx/myCGSEM % h(i,j,iEl)**2 - myCGSEM % f(i,j,iEl)/myCGSEM % h(i,j,iEl)*dudy
               !Fy = myCGSEM % params % cDrag*dudy/myCGSEM % h(i,j,iEl)**2 + myCGSEM % f(i,j,iEl)/myCGSEM % h(i,j,iEl)*dudx
               
               Fx = myCGSEM % params % cDrag*dudx - myCGSEM % f(i,j,iEl)/myCGSEM % h(i,j,iEl)*dudy
               Fy = myCGSEM % params % cDrag*dudy + myCGSEM % f(i,j,iEl)/myCGSEM % h(i,j,iEl)*dudx
            
               
               F1(i,j) = ( Fx*myCGSEM % mesh % geom(iEl) % dydp(i,j) - &
                           Fy*myCGSEM % mesh % geom(iEl) % dxdp(i,j) )*&
                           myCGSEM % cgStorage % qweight(i)
                           
               F2(i,j) = ( Fy*myCGSEM % mesh % geom(iEl) % dxds(i,j) - &
                           Fx*myCGSEM % mesh % geom(iEl) % dyds(i,j) )*&
                           myCGSEM % cgStorage % qweight(j)        

            ENDDO !
         ENDDO 

         ! Flux divergence
         DO j = 0, myCGSEM % N
            DO i = 0, myCGSEM % N
            
               dF1ds = 0.0_prec
               dF2dp = 0.0_prec
               DO k = 0, myCGSEM % N
                  dF1ds = dF1ds + myCGSEM % cgStorage % dMatP(i,k)*F1(k,j) !MATMUL( myCGSEM % cgStorage % dMatP, F1 )
                  dF2dp = dF2dp + F2(i,k)*myCGSEM % cgStorage % dMatS(k,j) !MATMUL( F2, myCGSEM % cgStorage % dMatS )
               ENDDO
               
               Au(i,j,iEl) = -( dF1ds*myCGSEM % cgStorage % qweight(j) + dF2dp*myCGSEM % cgStorage % qweight(i) )
            ENDDO 
         ENDDO 
      
      ENDDO
!$OMP END DO

!$OMP END PARALLEL

 END SUBROUTINE MatrixAction_GeoBarS
!
!
!
 SUBROUTINE Residual_GeoBarS( myCGSEM, r )
 ! S/R Residual
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(inout) :: myCGSEM
   REAL(prec), INTENT(out)         :: r(0:myCGSEM % N, &
                                        0:myCGSEM % N, &
                                        1:myCGSEM % mesh % nElems )
   ! LOCAL
   INTEGER    :: iEl, i, j
   REAL(prec) :: Au(0:myCGSEM % N, &
                    0:myCGSEM % N, &
                    1:myCGSEM % mesh % nElems )
                    
                    
      CALL myCGSEM % MatrixAction( myCGSEM % psi, Au )
      !$OMP PARALLEL
      !$OMP DO
      DO iEl = 1, myCGSEM % mesh % nElems
         DO j = 0, myCGSEM % N
            DO i = 0, myCGSEM % N
              
               r(i,j,iEl) = myCGSEM % source(i,j,iEl)*myCGSEM % mesh % geom(iEl) % J(i,j)*&
                            myCGSEM % cgStorage % qWeight(i)*myCGSEM % cgStorage % qWeight(j) -&
                            Au(i,j,iEl)
            ENDDO
         ENDDO
      ENDDO
      !$OMP ENDDO
      !$OMP END PARALLEL

     ! Add the contributions from the shared edges and corners
     CALL myCGSEM % GlobalSum( r )
     CALL myCGSEM % Mask( r )

 END SUBROUTINE Residual_GeoBarS
!
!
!
 FUNCTION DotProduct_GeoBarS( myCGSEM, u, v ) RESULT( uDotv )
 ! FUNCTION DotProduct
 ! Computes the dot product of two "vectors" whose data are
 ! stored in the "(iS,iP)" format
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( GeoBarS ) :: myCGSEM
   REAL(prec)       :: u(0:myCGSEM % N, &
                                  0:myCGSEM % N, &
                                  1:myCGSEM % mesh % nElems )
   REAL(prec)        :: v(0:myCGSEM % N, &
                                  0:myCGSEM % N, &
                                  1:myCGSEM % mesh % nElems )
   REAL(prec)        :: uDotv
   ! Local
   INTEGER :: iS, iP, N, iEl, nEl, k
   REAL(prec) :: uloc(0:myCGSEM % N, &
                      0:myCGSEM % N, &
                      1:myCGSEM % mesh % nElems )
   REAL(prec) :: vloc(0:myCGSEM % N, &
                      0:myCGSEM % N, &
                      1:myCGSEM % mesh % nElems )
   
      N   = myCGSEM % N
      nEl = myCGSEM % mesh % nElems
      
      uloc = u
      vloc = v
      ! mask the arrays
      CALL myCGSEM % Mask( uloc )
      CALL myCGSEM % Mask( vloc )
 
      uDotv = ZERO
      DO iEl = 1, nEl
         DO iP = 0, N
            DO iS = 0, N
               uDotV = uDotV + uloc(iS,iP,iEl)*vloc(iS,iP,iEl)
            ENDDO
         ENDDO
      ENDDO
      

 END FUNCTION DotProduct_GeoBarS
!
!
!
 SUBROUTINE Solve_GMRES( myCGSEM, ioerr )
 !  S/R Solve
 !
 !  This subroutine solves the system Ax = b using the un-preconditioned GMRES.
 !  The matrix action and residual routines are supplied by a non-abstracted type-extension of
 !  GMRES. These routines should return an array indexed from 1 to nDOF. Thus,
 !  in addition to a MatrixAction and Residual, the user should map their data-structure to a 1-D 
 !  array.  
 !
 !  On output ioerr is set to an error checking flag. 
 !  If ioerr ==  0, the method converged within the maximum number of iterations.
 !     ioerr == -1, the method did not converge within the maximum number of iterations.
 !     ioerr == -2, something that is not caught by the current construct happened.
 !  
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(inout) :: myCGSEM
   INTEGER, INTENT(out)            :: ioerr
   ! LOCAL
   INTEGER    :: i, j, k, l, iS, iP, iEl
   INTEGER    :: nIt, m, nr 
   REAL(prec) :: TOL
   REAL(prec) :: r(0:myCGSEM % N, 0:myCGSEM % N, 1:myCGSEM % mesh % nElems)
   REAL(prec) :: v(0:myCGSEM % N, 0:myCGSEM % N, 1:myCGSEM % mesh % nElems, 1:myCGSEM % params % mInnerIters+1)
   REAL(prec) :: w(0:myCGSEM % N, 0:myCGSEM % N, 1:myCGSEM % mesh % nElems)
   REAL(prec) :: w0(0:myCGSEM % N, 0:myCGSEM % N, 1:myCGSEM % mesh % nElems)
   REAL(prec) :: rho(1:myCGSEM % params % mInnerIters, 1:myCGSEM % params % mInnerIters)
   REAL(prec) :: h(1:myCGSEM % params % mInnerIters+1, 1:myCGSEM % params % mInnerIters)
   REAL(prec) :: c(1:myCGSEM % params % mInnerIters), s(1:myCGSEM % params % mInnerIters), bhat(1:myCGSEM % params % mInnerIters+1)
   REAL(prec) :: y(1:myCGSEM % params % mInnerIters+1)
   REAL(prec) :: b, d, g, r0, rc
   
      ioerr = -2
      nIt = myCGSEM % params % maximumIterates
      m   = myCGSEM % params % mInnerIters
      TOL = myCGSEM % tolerance
      r = ZERO

      CALL myCGSEM % Residual( r )
      
      r0 = sqrt(myCGSEM % DotProduct( r, r ) )
     
      l = 0 
      myCGSEM % resi = ZERO
      myCGSEM % resi(l) = r0
      
      v    = ZERO  
      rho  = ZERO
      bhat = ZERO 
      s    = ZERO
      c    = ZERO
      y    = ZERO
      
     
      DO j = 1,nIt

         b       = sqrt( myCGSEM % DotProduct( r, r ) )
         v(:,:,:,1)  = r/b
         bhat(1) = b

         DO i = 1, m
            l = l+1
            nr = i
            ! The first step in GMRES is to build the orthogonal basis up to order "i"
            ! with the accompanying upper hessenburg matrix.
            CALL myCGSEM % UnMask( v(:,:,:,i) )
            CALL myCGSEM % MatrixAction( v(:,:,:,i), w )
            CALL myCGSEM % GlobalSum( w )
            CALL myCGSEM % Mask( w )
            CALL myCGSEM % Mask( v(:,:,:,i) )
            
            ! The new basis vector is obtained by multiplying the previous basis vector by the matrix
            ! and orthogonalizing wrt to all of the previous basis vectors using a Gram-Schmidt process.
            DO k = 1, i
               h(k,i) = myCGSEM % DotProduct( v(:,:,:,k), w )
               w      = w - h(k,i)*v(:,:,:,k)
            ENDDO

            h(i+1,i) = sqrt( myCGSEM % DotProduct(w,w) )

            IF( AlmostEqual( h(i+1,i), ZERO )  )THEN
               EXIT
            ENDIF

            v(:,:,:,i+1) = w/h(i+1,i)
            rho(1,i) = h(1,i)

            ! Givens rotations are applied to the upper hessenburg matrix and to the residual vectors
            ! that are formed from the orthogonalization process. Here, they are done "on-the-fly"
            ! as opposed to building the entire upper hessenburg matrix and orthonormal basis
            ! before performing the rotations. This way, we can also tell if we have found an exact
            ! solution ( if h(i+1,i) = 0 ) with a smaller subspace than size m.
            DO k = 2, i

               g          = c(k-1)*rho(k-1,i) + s(k-1)*h(k,i)
               rho(k,i)   = -s(k-1)*rho(k-1,i) + c(k-1)*h(k,i)
               rho(k-1,i) = g 

            ENDDO

            ! Here the coefficients of the Givens rotation matrix are computed
            d = sqrt( rho(i,i)**2 + h(i+1,i)**2 )
            c(i) = rho(i,i)/d
            s(i) = h(i+1,i)/d

            rho(i,i) = c(i)*rho(i,i) + s(i)*h(i+1,i)
            ! And applied to the residual vector
            bhat(i+1) = -s(i)*bhat(i)
            bhat(i)   = c(i)*bhat(i)
 
         
            rc = abs( bhat(i+1) )

            myCGSEM % resi(l) = rc
    !        print*, rc
            IF( rc <= TOL )THEN
               EXIT
            ENDIF

         ENDDO
         
         IF( rc > TOL )THEN
            nr = m
         ENDIF

         ! Back-substitution of the tridiagonal matrix that resulted from the rotations
         y(nr) = bhat(nr)/rho(nr,nr)
         DO k = nr-1, 1, -1

            y(k) = bhat(k)

            DO i = k+1, nr
               y(k) = y(k) - rho(k,i)*y(i)

            ENDDO

            y(k) = y(k)/rho(k,k)

         ENDDO
         
 
         !myCGSEM % psi = ZERO
         DO iEl = 1,myCGSEM % mesh % nElems
            DO iP = 0, myCGSEM % N
               DO iS = 0, myCGSEM % N
                  myCGSEM % psi(iS,iP,iEl) = myCGSEM % psi(iS,iP,iEl) + &
                                             DOT_PRODUCT( v(iS,iP,iEl,1:nr), y(1:nr) )
               ENDDO
            ENDDO
         ENDDO
         CALL myCGSEM % UnMask( myCGSEM % psi )
         
         IF( rc <= TOL )THEN
            ioerr = 0
            EXIT
         ENDIF
         CALL myCGSEM % Residual( r )
         

      ENDDO 

      IF( rc > TOL )THEN
         PRINT*, 'MODULE IterativeSolvers : GMRES failed to converge '
         PRINT*, 'Last L-2 residual : ', sqrt(myCGSEM % DotProduct(r,r))
         ioerr=-1
      ENDIF   

 END SUBROUTINE Solve_GMRES
!
!
!
SUBROUTINE DiagnoseVelocity_GeoBarS( myCGSEM )
 ! S/R DiagnoseVelocity
 ! 
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS(GeoBarS), INTENT(inout) :: myCGSEM
   ! Local
   INTEGER :: i, j, k, iEl
   REAL(prec) :: dpds, dpdp, dpdx, dpdy


     
!$OMP PARALLEL


!$OMP DO PRIVATE( dpds, dpdp, dpdx, dpdy )
      DO iEl = 1, myCGSEM % mesh % nElems
      
      
         !------------- Contravariant Flux calculations ---------------!
         ! Calculate the gradient in the stream function 
         !CALL myCGSEM % CalculateGradient( iEl, u, dudx, dudy ) ! calculate the gradient in physical coordinates
         DO j = 0, myCGSEM % N
            DO i = 0, myCGSEM % N
            
               dpds = 0.0_prec
               dpdp = 0.0_prec
               
               ! Gradient in computational space
               DO k = 0, myCGSEM % N
                  dpds = dpds + myCGSEM % cgStorage % dMatS(i,k)*myCGSEM % psi(k,j,iEl)
                  dpdp = dpdp + myCGSEM % psi(i,k,iEl)*myCGSEM % cgStorage % dMatP(k,j)
               ENDDO
               
               
               ! Calculate the gradient in physical space
               dpdx = ( myCGSEM % mesh % geom(iEl) % dydp(i,j)*dpds - &
                        myCGSEM % mesh % geom(iEl) % dyds(i,j)*dpdp )/&
                        myCGSEM % mesh % geom(iEl) % J(i,j)
                        
               dpdy = ( myCGSEM % mesh % geom(iEl) % dxds(i,j)*dpdp - &
                        myCGSEM % mesh % geom(iEl) % dxdp(i,j)*dpds )/&
                        myCGSEM % mesh % geom(iEl) % J(i,j)

               myCGSEM % u(i,j,iEl) =  dpdy
               myCGSEM % v(i,j,iEl) = -dpdx
               
            ENDDO !
         ENDDO 

      
      ENDDO
!$OMP END DO

!$OMP END PARALLEL

 END SUBROUTINE DiagnoseVelocity_GeoBarS

!
!
!==================================================================================================!
!---------------------------------------- File I/O ------------------------------------------------!
!==================================================================================================!
!
!
SUBROUTINE CoarseToFine_GeoBarS( myCGSEM, iEl, x, y, sol, f, h, u, v, source )
 ! CoarseToFine
 !
 ! =============================================================================================== !
 ! DECLARATIONS
   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(in) :: myCGSEM
   INTEGER, INTENT(in)             :: iEl
   REAL(prec), INTENT(out)         :: x(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)         :: y(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)         :: sol(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)         :: f(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)         :: h(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)         :: u(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)         :: v(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec), INTENT(out)         :: source(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   ! Local
   INTEGER    :: i, nModes
   
      
      x   = myCGSEM % cgStorage % interp % ApplyInterpolationMatrix_2D( myCGSEM % mesh % geom(iEl) % x )
      y   = myCGSEM % cgStorage % interp % ApplyInterpolationMatrix_2D( myCGSEM % mesh % geom(iEl) % y )
      sol = myCGSEM % cgStorage % interp % ApplyInterpolationMatrix_2D( myCGSEM % psi(:,:,iEl) )
      f   = myCGSEM % cgStorage % interp % ApplyInterpolationMatrix_2D( myCGSEM % f(:,:,iEl) )
      h   = myCGSEM % cgStorage % interp % ApplyInterpolationMatrix_2D( myCGSEM % h(:,:,iEl) )
      u   = myCGSEM % cgStorage % interp % ApplyInterpolationMatrix_2D( myCGSEM % u(:,:,iEl) )
      v   = myCGSEM % cgStorage % interp % ApplyInterpolationMatrix_2D( myCGSEM % v(:,:,iEl) )
      source = myCGSEM % cgStorage % interp % ApplyInterpolationMatrix_2D( myCGSEM % source(:,:,iEl) )
      
 END SUBROUTINE CoarseToFine_GeoBarS
!
!
!
 SUBROUTINE WriteTecplot_GeoBarS( myCGSEM, filename )
 ! S/R WriteTecplot
 !
 !
 ! =============================================================================================== !
   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(inout) :: myCGSEM
   CHARACTER(*), INTENT(in), OPTIONAL :: filename
   !LOCAL
   INTEGER :: i, j, iEl, fUnit
   CHARACTER(len=5) :: zoneID
   REAL(prec)       :: x(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: y(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: sol(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: h(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: f(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: u(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: v(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)
   REAL(prec)       :: source(0:myCGSEM % nPlot, 0:myCGSEM % nPlot)

         OPEN( UNIT=NewUnit(fUnit), &
               FILE='GeoBarS.tec', &
               FORM='FORMATTED')
      
      WRITE(fUnit,*) 'VARIABLES = "X", "Y","psi","f","h","u","v" "source"'

      CALL myCGSEM % DiagnoseVelocity( )
      
      DO iEl = 1, myCGSEM % mesh % nElems

         WRITE(zoneID,'(I5.5)') iEl
         WRITE(fUnit,*) 'ZONE T="el'//trim(zoneID)//'", I=',myCGSEM % params % nPlot+1,', J=', myCGSEM % params % nPlot+1,',F=POINT'
         CALL myCGSEM % CoarseToFine( iEl, x, y, sol, f, h, u, v, source )

         DO j = 0, myCGSEM % params % nPlot
            DO i = 0, myCGSEM % params % nPlot
               WRITE(fUnit,*)  x(i,j), y(i,j), sol(i,j), f(i,j), h(i,j), u(i,j), v(i,j), source(i,j)
            ENDDO
         ENDDO

      ENDDO
      
      CLOSE( fUnit )
      

 END SUBROUTINE WriteTecplot_GeoBarS
!
!
!
 SUBROUTINE WriteResidual_GeoBarS( myCGSEM )
 
   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(in) :: myCGSEM
   !  Local
   INTEGER :: i, fUnit
   
   
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE='Residual.curve', &
            FORM='formatted',&
            ACCESS='sequential',&
            STATUS='replace',&
            ACTION='WRITE')
            
            
       DO i = 1, myCGSEM % nresi
          WRITE(fUnit,*) i, myCGSEM % resi(i)
       ENDDO
       
       CLOSE(fUnit)
   
 END SUBROUTINE WriteResidual_GeoBarS
!
!
!
 SUBROUTINE WritePickup_GeoBarS(myCGSEM)

   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(in) :: myCGSEM
  ! LOCAL
   INTEGER       :: i
   INTEGER       :: thisRec, fUnit
   INTEGER       :: N, nElems

      N      = myCGSEM % N
      nElems = myCGSEM % mesh % nElems
      
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE='GeoBarS.pickup', &
            FORM='unformatted',&
            ACCESS='direct',&
            STATUS='replace',&
            ACTION='WRITE',&
            CONVERT='big_endian',&
            RECL=prec*(N+1)*(N+1)*nElems )

      thisRec = 1 
      WRITE( fUnit, REC=thisRec ) myCGSEM % psi(0:N,0:N,1:nElems) 
      thisRec = thisRec+1
      WRITE( fUnit, REC=thisRec ) myCGSEM % h(0:N,0:N,1:nElems) 
      thisRec = thisRec+1
      WRITE( fUnit, REC=thisRec ) myCGSEM % f(0:N,0:N,1:nElems) 
      thisRec = thisRec+1

      CLOSE(UNIT=fUnit)


 END SUBROUTINE WritePickup_GeoBarS
!
!
 SUBROUTINE ReadPickup_GeoBarS(myCGSEM)

   IMPLICIT NONE
   CLASS( GeoBarS ), INTENT(inout) :: myCGSEM
  ! LOCAL
   INTEGER       :: i
   INTEGER       :: thisRec, fUnit
   INTEGER       :: N, nElems

      N      = myCGSEM % N
      nElems = myCGSEM % mesh % nElems
      
      OPEN( UNIT=NEWUNIT(fUnit), &
            FILE='GeoBarS.pickup', &
            FORM='unformatted',&
            ACCESS='direct',&
            STATUS='OLD',&
            ACTION='READ',&
            CONVERT='big_endian',&
            RECL=prec*(N+1)*(N+1)*nElems )

      thisRec = 1 
      READ( fUnit, REC=thisRec ) myCGSEM % psi(0:N,0:N,1:nElems) 
      thisRec = thisRec+1
      READ( fUnit, REC=thisRec ) myCGSEM % h(0:N,0:N,1:nElems) 
      thisRec = thisRec+1
      READ( fUnit, REC=thisRec ) myCGSEM % f(0:N,0:N,1:nElems) 
      thisRec = thisRec+1
      
      
      CLOSE(UNIT=fUnit)


 END SUBROUTINE ReadPickup_GeoBarS
!
END MODULE GeoBarS_Class
