! ModelPrecision.f90
! 
! Copyright 2018 Joseph Schoonover <joeschoonover@fluidnumerics.com>, Fluid Numerics, LLC
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

!> \file ModelPrecision.f90
!! Contains the \ref ModelPrecision module

!> \defgroup ModelPrecision ModelPrecision 
!! This module is used to set the precision of floating point numbers throughout the rest of the code.
!!
!! To set the precision to <B>single precision</B>, set <B>prec=sp</B> <BR>
!! To set the precision to <B>double precision</B>, set <B>prec=dp</B> <BR>
!! The precision must be set prior to compilation.

 
MODULE ModelPrecision

INTEGER, PARAMETER :: sp   = SELECTED_REAL_KIND(6, 37)     ! 32-bit
INTEGER, PARAMETER :: dp   = SELECTED_REAL_KIND(15, 307)   ! 64-bit
INTEGER, PARAMETER :: qp   = SELECTED_REAL_KIND(33, 4931)  ! 128-bit
INTEGER, PARAMETER :: prec = dp                            ! Specify the precision here

END MODULE ModelPrecision
